#!/usr/bin/python3

import asyncio
import json
import logging
import os
import re
import tempfile
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple
from aiogram.types import BotCommand

import requests
import redis

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from aiogram import Bot, Dispatcher, F
from aiogram.enums import ChatAction
from aiogram.filters import Command, CommandObject
from aiogram.types import (
    CallbackQuery,
    InlineKeyboardButton,
    InlineKeyboardMarkup,
    Message,
    FSInputFile,
    BotCommand
)

from dotenv import load_dotenv

load_dotenv()

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

async def set_main_menu(bot: Bot):

    commands = [
        BotCommand(command="/start", description="Начало работы"),
        BotCommand(command="/help", description="Помощь и справка"),
        BotCommand(command="/get", description="Получить данные по rsID"),
        BotCommand(command="/history", description="История запросов"),
        BotCommand(command="/about", description="О боте")
    ]
    await bot.set_my_commands(commands)

@dataclass
class Config:
    TELEGRAM_BOT_TOKEN: str
    REDIS_URL: str
    NCBI_API_BASE: str
    NCBI_API_TIMEOUT: int
    CACHE_TTL: int
    MAX_REQUESTS_PER_HOUR: int
    EXAMPLE_RSIDS: Tuple[str, ...]

def load_config() -> Config:
    return Config(
        TELEGRAM_BOT_TOKEN=os.getenv("TELEGRAM_BOT_TOKEN", "").strip(),
        REDIS_URL=os.getenv("REDIS_URL", "redis://localhost:6379/0"),
        NCBI_API_BASE=os.getenv(
            "NCBI_API_BASE",
            "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp",
        ),
        NCBI_API_TIMEOUT=int(os.getenv("NCBI_API_TIMEOUT", "30")),
        CACHE_TTL=int(os.getenv("CACHE_TTL", "86400")),
        MAX_REQUESTS_PER_HOUR=int(os.getenv("MAX_REQUESTS_PER_HOUR", "50")),
        EXAMPLE_RSIDS=tuple(
            os.getenv("EXAMPLE_RSIDS", "rs7755898,rs429358,rs7412").split(",")
        ),
    )

class CacheManager:
    def __init__(self, config: Config):
        self.config = config
        self._redis: Optional[redis.Redis] = None
        self._init_redis()

    def _init_redis(self) -> None:
        try:
            self._redis = redis.Redis.from_url(
                self.config.REDIS_URL,
                decode_responses=True,
            )
            self._redis.ping()
            logger.info("Connected to Redis at %s", self.config.REDIS_URL)
        except Exception as e:
            logger.warning("Redis is not available: %s", e)
            self._redis = None

    def get_snp(self, rsid: str) -> Optional[Dict[str, Any]]:
        if not self._redis:
            return None
        key = f"snp:{rsid.lower()}"
        raw = self._redis.get(key)
        if not raw:
            return None
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            logger.warning("Failed to decode cached SNP JSON for %s", rsid)
            return None

    def set_snp(self, rsid: str, data: Dict[str, Any]) -> None:
        if not self._redis:
            return
        key = f"snp:{rsid.lower()}"
        self._redis.setex(key, self.config.CACHE_TTL, json.dumps(data))

    def add_history(self, user_id: int, rsid: str, max_items: int = 10) -> None:
        if not self._redis:
            return
        key = f"history:{user_id}"
        pipe = self._redis.pipeline()
        pipe.lpush(key, rsid)
        pipe.ltrim(key, 0, max_items - 1)
        pipe.expire(key, self.config.CACHE_TTL)
        pipe.execute()

    def get_history(self, user_id: int, limit: int = 10) -> List[str]:
        if not self._redis:
            return []
        key = f"history:{user_id}"
        return self._redis.lrange(key, 0, limit - 1)

    def check_rate_limit(self, user_id: int) -> Tuple[bool, Optional[int]]:
        limit = self.config.MAX_REQUESTS_PER_HOUR
        if limit <= 0 or not self._redis:
            return True, None

        now = datetime.now(timezone.utc)
        key_hour = now.strftime("%Y%m%d%H")
        key = f"rate:{user_id}:{key_hour}"

        try:
            pipe = self._redis.pipeline()
            pipe.incr(key)
            pipe.expire(key, 3600)
            count, _ = pipe.execute()
            allowed = count <= limit
            return allowed, int(count)
        except Exception as e:
            logger.warning("Rate limiting failed: %s", e)
            return True, None

@dataclass
class AlleleFrequency:
    allele: str
    study: str
    frequency: float
    allele_count: int
    total_count: int

@dataclass
class GenotypeFrequencies:
    study: str
    ref_allele: str
    alt_allele: str
    p_ref: float
    q_alt: float
    p2_hom_ref: float
    two_pq_het: float
    q2_hom_alt: float
    total_alleles: int

@dataclass
class SNPAnalysisResult:
    rsid: str
    ref_snp_id: int
    allele_frequencies: List[AlleleFrequency]
    genotype_frequencies: List[GenotypeFrequencies]
    raw_json: Dict[str, Any]

def fetch_snp_json(config: Config, rsid: str) -> Dict[str, Any]:
    if not rsid.lower().startswith("rs"):
        raise ValueError("rsID должен начинаться с 'rs'")

    try:
        numeric_id = int(rsid[2:])
    except ValueError:
        raise ValueError("После 'rs' должны быть только цифры")

    url = f"{config.NCBI_API_BASE}/{numeric_id}"
    logger.info("Requesting NCBI API: %s", url)
    resp = requests.get(url, timeout=config.NCBI_API_TIMEOUT)

    if resp.status_code == 404:
        raise ValueError(f"rsID {rsid} не найден в базе данных NCBI dbSNP.")

    if not resp.ok:
        raise RuntimeError(
            f"Ошибка при обращении к NCBI API: {resp.status_code} {resp.text[:200]}"
        )

    return resp.json()

def extract_allele_frequencies(data: Dict[str, Any]) -> List[AlleleFrequency]:
    result: List[AlleleFrequency] = []

    primary = data.get("primary_snapshot_data", {})
    allele_annotations = primary.get("allele_annotations", [])

    for allele_ann in allele_annotations:
        freq_list = allele_ann.get("frequency", [])
        for freq_record in freq_list:
            try:
                obs = freq_record.get("observation", {})
                deleted_seq = obs.get("deleted_sequence", "")
                inserted_seq = obs.get("inserted_sequence", "")
                allele_count = int(freq_record.get("allele_count", 0))
                total_count = int(freq_record.get("total_count", 0))
                study_name = freq_record.get("study_name", "Unknown")

                if total_count <= 0 or allele_count < 0:
                    continue

                freq_val = allele_count / total_count
                allele = inserted_seq if inserted_seq else deleted_seq
                if not allele:
                    continue

                result.append(
                    AlleleFrequency(
                        allele=allele,
                        study=study_name,
                        frequency=freq_val,
                        allele_count=allele_count,
                        total_count=total_count,
                    )
                )
            except Exception as e:
                logger.warning("Failed to parse frequency record: %s", e)

    return result

def calculate_hardy_weinberg(
    alleles: List[AlleleFrequency],
) -> List[GenotypeFrequencies]:
    by_study: Dict[str, Dict[str, AlleleFrequency]] = {}

    for af in alleles:
        by_study.setdefault(af.study, {})
        existing = by_study[af.study].get(af.allele)
        if existing is None or af.frequency > existing.frequency:
            by_study[af.study][af.allele] = af

    results: List[GenotypeFrequencies] = []

    for study, allele_map in by_study.items():
        if len(allele_map) < 2:
            continue

        sorted_alleles = sorted(
            allele_map.values(), key=lambda x: x.frequency, reverse=True
        )
        a1, a2 = sorted_alleles[0], sorted_alleles[1]

        total = a1.frequency + a2.frequency
        if total == 0:
            continue

        p = a1.frequency / total
        q = a2.frequency / total

        p2 = p * p
        two_pq = 2 * p * q
        q2 = q * q

        total_alleles = max(a1.total_count, a2.total_count)

        results.append(
            GenotypeFrequencies(
                study=study,
                ref_allele=a1.allele,
                alt_allele=a2.allele,
                p_ref=p,
                q_alt=q,
                p2_hom_ref=p2,
                two_pq_het=two_pq,
                q2_hom_alt=q2,
                total_alleles=total_alleles,
            )
        )

    return results

def format_text_report(result: SNPAnalysisResult) -> str:
    lines: List[str] = []

    lines.append(f"Результаты для запроса {result.rsid}")
    lines.append("")
    lines.append(f"Референсный SNP ID: {result.ref_snp_id}")
    lines.append("")

    lines.append("Найденные варианты аллелей:")
    lines.append("")

    by_allele: Dict[str, List[AlleleFrequency]] = {}
    for af in result.allele_frequencies:
        by_allele.setdefault(af.allele, []).append(af)

    for allele, freq_list in by_allele.items():
        lines.append(f"🔸 Аллель: {allele}")
        for af in freq_list:
            lines.append(
                f"  • {af.study}: {af.frequency:.6f} "
                f"({af.allele_count}/{af.total_count})"
            )
        lines.append("")

    lines.append("Рассчитанные частоты генотипов:")
    lines.append("")

    for gf in result.genotype_frequencies:
        lines.append(f"🔸 Исследование: {gf.study}")
        lines.append(
            f"  Референсный аллель (ref, 0): {gf.ref_allele} (частота: {gf.p_ref:.6f})"
        )
        lines.append(
            f"  Альтернативный аллель (alt, 1): {gf.alt_allele} (частота: {gf.q_alt:.6f})"
        )
        lines.append("  Частоты генотипов:")
        lines.append(f"    ├── 0/0 (гомозиготный REF): {gf.p2_hom_ref:.6f}")
        lines.append(f"    ├── 0/1 (гетерозиготный): {gf.two_pq_het:.6f}")
        lines.append(f"    └── 1/1 (гомозиготный ALT): {gf.q2_hom_alt:.6f}")
        lines.append(f"  Всего аллелей: {gf.total_alleles}")
        lines.append("")
        lines.append("")

    lines.append(
        f"Всего исследований: {len(result.genotype_frequencies)}"
    )

    return "\n".join(lines)

def analyze_snp(config: Config, rsid: str) -> SNPAnalysisResult:
    data = fetch_snp_json(config, rsid)

    ref_snp_id = data.get("refsnp_id", 0)
    try:
        ref_snp_id = int(ref_snp_id)
    except Exception:
        ref_snp_id = 0

    allele_freqs = extract_allele_frequencies(data)
    if not allele_freqs:
        raise RuntimeError("Не удалось извлечь частоты аллелей из ответа NCBI.")

    geno_freqs = calculate_hardy_weinberg(allele_freqs)

    return SNPAnalysisResult(
        rsid=rsid,
        ref_snp_id=ref_snp_id,
        allele_frequencies=allele_freqs,
        genotype_frequencies=geno_freqs,
        raw_json=data,
    )

def to_json_serializable(result: SNPAnalysisResult) -> Dict[str, Any]:
    return {
        "rsid": result.rsid,
        "ref_snp_id": result.ref_snp_id,
        "allele_frequencies": [asdict(a) for a in result.allele_frequencies],
        "genotype_frequencies": [asdict(g) for g in result.genotype_frequencies],
        "raw_json": result.raw_json,
    }

def dump_json(result: SNPAnalysisResult, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(to_json_serializable(result), f, ensure_ascii=False, indent=2)

def _choose_study_for_plot(result: SNPAnalysisResult) -> str:
    names = [g.study for g in result.genotype_frequencies]
    if not names:
        raise RuntimeError("Нет генотипических частот для построения графиков.")

    for preferred in ("1000Genomes_30X", "1000Genomes", "GnomAD_exomes", "GnomAD_genomes"):
        if preferred in names:
            return preferred

    return names[0]

def generate_plots(
    result: SNPAnalysisResult,
    output_dir: str,
) -> Tuple[str, str]:
    os.makedirs(output_dir, exist_ok=True)
    study_name = _choose_study_for_plot(result)

    allele_freqs: Dict[str, float] = {}
    total_alleles = 0
    for af in result.allele_frequencies:
        if af.study != study_name:
            continue
        allele_freqs[af.allele] = af.frequency
        total_alleles = max(total_alleles, af.total_count)

    if not allele_freqs:
        raise RuntimeError(
            f"Нет аллельных частот для исследования {study_name}."
        )
    
    plt.style.use('bmh')
    
    alleles = list(allele_freqs.keys())
    freqs = [allele_freqs[a] for a in alleles]
    
    fig, ax = plt.subplots()
    ax.bar(alleles, freqs, color = 'yellowgreen')
    
    ax.set_title(f"Распределение аллельных частот для ({study_name})", fontweight='bold', fontsize=11, color = 'darkgreen')
    ax.set_xlabel("Аллель",  fontsize=11, color = 'darkgreen')
    ax.set_ylabel("Частота",  fontsize=11, color = 'darkgreen')
    
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    
    for i, freq in enumerate(freqs):
        ax.text(i, freq + 0.01, f"{freq:.3f}", 
                ha='center', va='bottom', fontsize=8)
    
    allele_plot_path = os.path.join(
        output_dir, f"{result.rsid}_{study_name}_alleles.png"
    )
    
    fig.tight_layout()
    fig.savefig(allele_plot_path, dpi=450)
    plt.close(fig)

    gf = next(g for g in result.genotype_frequencies if g.study == study_name)
    
    labels = [
        f"{gf.ref_allele}/{gf.ref_allele}",
        f"{gf.ref_allele}/{gf.alt_allele}",
        f"{gf.alt_allele}/{gf.alt_allele}",
    ]
    sizes = [gf.p2_hom_ref, gf.two_pq_het, gf.q2_hom_alt]
    
    fig2, ax2 = plt.subplots()
    explode = (0.05, 0.05, 0.05) 
    
    wedges, texts, autotexts = ax2.pie(
        sizes,
        labels=labels,
        autopct="%1.2f",
        explode=explode,
        startangle=90,
        shadow=False,
        textprops={ 'fontsize': 5}
    )
    
    for autotext in autotexts:
        autotext.set_fontsize(5)
    
    ax2.set_title(f"Распределение генотипов у ({study_name})", fontweight='bold', color = 'darkred', fontsize=11)
    
    ax2.legend(
        wedges,
        labels,
        title="Генотипы",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1),
        prop={'size': 9}
    )
    
    ax2.axis('equal')
    
    geno_plot_path = os.path.join(
        output_dir, f"{result.rsid}_{study_name}_genotypes.png"
    )
    
    fig2.tight_layout()
    fig2.savefig(geno_plot_path, dpi=450)
    plt.close(fig2)
    
    logger.info("Plots generated: %s, %s", allele_plot_path, geno_plot_path)
    return allele_plot_path, geno_plot_path

RSID_REGEX = re.compile(r"^rs\d+$", re.IGNORECASE)

def validate_rsid(rsid: str) -> bool:
    return bool(RSID_REGEX.match(rsid.strip()))

def build_example_keyboard(example_rsids: Tuple[str, ...]) -> InlineKeyboardMarkup:
    buttons = [
        [
            InlineKeyboardButton(
                text=rsid.strip(),
                callback_data=f"example:{rsid.strip()}",
            )
        ]
        for rsid in example_rsids
        if rsid.strip()
    ]
    return InlineKeyboardMarkup(inline_keyboard=buttons)

async def process_rsid(
    bot: Bot,
    config: Config,
    cache: CacheManager,
    rsid: str,
    chat_id: int,
    user_id: int,
    progress_message: Optional[Message] = None,
) -> None:
    rsid = rsid.strip()

    if not validate_rsid(rsid):
        await bot.send_message(
            chat_id,
            "Неверный формат rsID. Используйте вид <code>rs123456</code>.",
            parse_mode="HTML"
        )
        return

    allowed, count = cache.check_rate_limit(user_id)
    if not allowed:
        await bot.send_message(
            chat_id,
            "⚠️ <b>Превышен лимит запросов в час.</b> Попробуйте позже.",
            parse_mode="HTML"
        )
        return

    if progress_message is None:
        progress_message = await bot.send_message(
            chat_id, f"🔄 Получаю данные для {rsid}..."
        )
    else:
        await progress_message.edit_text(f"🔄 Получаю данные для {rsid}...")

    await bot.send_chat_action(chat_id, ChatAction.TYPING)

    cached = cache.get_snp(rsid)
    if cached:
        logger.info("Using cached result for %s", rsid)
        allele_freqs = [
            AlleleFrequency(**a) for a in cached.get("allele_frequencies", [])
        ]
        geno_freqs = [
            GenotypeFrequencies(**g) for g in cached.get("genotype_frequencies", [])
        ]
        result = SNPAnalysisResult(
            rsid=cached.get("rsid", rsid),
            ref_snp_id=cached.get("ref_snp_id", 0),
            allele_frequencies=allele_freqs,
            genotype_frequencies=geno_freqs,
            raw_json=cached.get("raw_json", {}),
        )
    else:
        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(None, analyze_snp, config, rsid)
        cache.set_snp(rsid, to_json_serializable(result))

    cache.add_history(user_id, rsid)

    with tempfile.TemporaryDirectory() as tmpdir:
        text_report = format_text_report(result)

        json_path = os.path.join(tmpdir, f"{rsid}_results.json")
        dump_json(result, json_path)

        allele_plot_path, geno_plot_path = generate_plots(result, tmpdir)

        await progress_message.edit_text(
            f"✅ Результаты для {rsid} готовы. Отправляю..."
        )

        if len(text_report) > 3800:
            text_for_tg = text_report[:3800] + "\n\n[Обрезано...]"
        else:
            text_for_tg = text_report

        await bot.send_message(chat_id, text_for_tg)

        await bot.send_photo(
            chat_id,
            FSInputFile(allele_plot_path),
            caption="📊 <b>Столбчатая диаграмма частот аллелей</b>",
            parse_mode="HTML"
        )
        await bot.send_photo(
            chat_id,
            FSInputFile(geno_plot_path),
            caption="📈 <b>Круговая диаграмма распределения генотипов</b>",
            parse_mode="HTML"
        )

        await bot.send_document(
            chat_id,
            FSInputFile(json_path),
            caption="📄 <b>Структурированные данные (JSON)</b>",
            parse_mode="HTML"
        )

async def main() -> None:
    config = load_config()
    if not config.TELEGRAM_BOT_TOKEN:
        raise RuntimeError("TELEGRAM_BOT_TOKEN не задан.")

    bot = Bot(token=config.TELEGRAM_BOT_TOKEN)
    dp = Dispatcher()
    cache = CacheManager(config)
    await set_main_menu(bot)

    @dp.message(Command("start"))
    async def cmd_start(message: Message) -> None:
        kb = build_example_keyboard(config.EXAMPLE_RSIDS)
        text = (
            "👋 <b>Привет! Я Sniphub — бот для анализа SNP</b>\n\n"
            
            "📊 <b>Основные команды:</b>\n"
            "• <code>/get rsID</code> — полный анализ SNP (частоты, графики)\n"
            "• <code>/history</code> — история запросов\n"
            "• <code>/help</code> — справка по командам\n"
            "• <code>/about</code> — информация о боте\n\n"
            
            "Используйте команду /get, например:\n"
            "  <code>/get rs7755898</code>\n\n"
            
            "👇 <i>Или нажми на одну из кнопок с примерами ниже.</i>"
        )
        await message.answer(text, reply_markup=kb, parse_mode="HTML")

    @dp.message(Command("help"))
    async def cmd_help(message: Message) -> None:
        text = (
            "<b>📚 Справка по Sniphub</b>\n\n"
            "🏠 /start — приветствие и краткая инструкция\n"
            "❓ /help — справка\n"
            "🔍 /get rsid — получить данные по rsID\n"
            "ℹ️ /about — информация о боте\n"
            "📜 /history — последние запросы\n\n"
            "<i>Пример:</i> <code>/get rs7755898</code>"
        )
        await message.answer(text, parse_mode="HTML")

    @dp.message(Command("about"))
    async def cmd_about(message: Message) -> None:
        text = (
            "<b>🧬 Sniphub — бот для анализа SNP</b>\n\n"
            
            "<b>📊 Функционал:</b>\n"
            "• Анализ частот аллелей из dbSNP\n"
            "• Расчет частот генотипов по Харди-Вайнбергу\n"
            "• Визуализация результатов\n\n"
            
            "<b>Источник данных:</b>\n"
            "• dbSNP: https://www.ncbi.nlm.nih.gov/snp/\n\n"
            
            "В случае обнаружения багов обращайтесь в поддержку телеграмма: @etcetera_mi"
        )
        await message.answer(text, parse_mode="HTML")

    @dp.message(Command("get"))
    async def cmd_get(message: Message, command: CommandObject) -> None:
        if not command.args:
            await message.answer(
                "❓ <b>Укажите rsID</b>, например:\n<code>/get rs7755898</code>",
                parse_mode="HTML"
            )
            return

        rsid = command.args.strip()
        progress_msg = await message.answer("🔄 <b>Запрос принят, начинаю обработку...</b>", parse_mode="HTML")
        try:
            await process_rsid(
                bot=bot,
                config=config,
                cache=cache,
                rsid=rsid,
                chat_id=message.chat.id,
                user_id=message.from_user.id,
                progress_message=progress_msg,
            )
        except Exception as e:
            logger.exception("Failed to process rsid %s: %s", rsid, e)
            await progress_msg.edit_text(
                f"❌ <b>Произошла ошибка при обработке {rsid}.</b>\nСообщение: {e}",
                parse_mode="HTML"
            )

    @dp.message(Command("history"))
    async def cmd_history(message: Message) -> None:
        history = cache.get_history(message.from_user.id)
        if not history:
            await message.answer("📭 <b>История пуста.</b>", parse_mode="HTML")
            return
        ordered = list(reversed(history))
        lines = ["<b>📜 История запросов по rsID:</b>", ""]
        for idx, rsid in enumerate(ordered, start=1):
            lines.append(f"{idx}. {rsid}")
        await message.answer("\n".join(lines), parse_mode="HTML")

    @dp.callback_query(F.data.startswith("example:"))
    async def on_example(callback: CallbackQuery) -> None:
        rsid = callback.data.split(":", 1)[1]
        await callback.answer(f"Обрабатываю запрос {rsid}...")
        msg = await callback.message.answer(
            f"🔄 <b>Запрос по rsID {rsid} выполнен.</b> Получаю данные и строю графики...",
            parse_mode="HTML"
        )
        try:
            await process_rsid(
                bot=bot,
                config=config,
                cache=cache,
                rsid=rsid,
                chat_id=callback.message.chat.id,
                user_id=callback.from_user.id,
                progress_message=msg,
            )
        except Exception as e:
            logger.exception("Failed to process example %s: %s", rsid, e)
            await msg.edit_text(
                f"❌ <b>Ошибка при обработке примера {rsid}:</b> {e}",
                parse_mode="HTML"
            )

    @dp.message()
    async def on_other(message: Message) -> None:
        await message.answer(
            "🤖 <b>Я понимаю только команды.</b>\n"
            "Используй <code>/get rsID</code>",
            parse_mode="HTML"
        )

    logger.info("Starting bot polling...")
    await dp.start_polling(bot)

asyncio.run(main())
