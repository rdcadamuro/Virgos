#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_tools.py  —  Bloco 1.a do pipeline viral.

Verifica se todas as ferramentas necessárias estão instaladas e funcionando.
Para cada uma:
  1) tenta detectar (comando no PATH + teste de versão);
  2) se faltar e --install for passado, tenta instalar na ordem de prioridade
     conda(bioconda) -> pip -> git  (cada fonte é um "alias" de fallback);
  3) imprime um relatório final CLARO do que está OK e do que ainda falta.

Uso:
  python check_tools.py --check
  python check_tools.py --install
  python check_tools.py --install --download-dbs --db-dir ~/viral_dbs

Projetado para WSL2/Linux com mamba ou conda no PATH.
"""

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Callable, List, Optional

# ----------------------------------------------------------------------
# Cores no terminal (degradam para vazio se não for TTY)
# ----------------------------------------------------------------------
def _c(code: str) -> str:
    return code if sys.stdout.isatty() else ""

GREEN, RED, YELLOW, BLUE, BOLD, RESET = (
    _c("\033[92m"), _c("\033[91m"), _c("\033[93m"),
    _c("\033[94m"), _c("\033[1m"), _c("\033[0m"),
)
OK   = f"{GREEN}✔{RESET}"
FAIL = f"{RED}x{RESET}"
WARN = f"{YELLOW}!{RESET}"


def run(cmd: List[str], quiet: bool = False) -> subprocess.CompletedProcess:
    """Executa um comando, captura saída, nunca levanta exceção."""
    if not quiet:
        print(f"  {BLUE}$ {' '.join(cmd)}{RESET}")
    return subprocess.run(cmd, capture_output=True, text=True)


def have(binname: str) -> bool:
    return shutil.which(binname) is not None


def conda_frontend() -> Optional[str]:
    """Prioriza mamba (mais rápido); cai para conda."""
    for fe in ("mamba", "conda"):
        if have(fe):
            return fe
    return None


# ----------------------------------------------------------------------
# Definição de cada ferramenta
# ----------------------------------------------------------------------
@dataclass
class Tool:
    name: str
    # comando(s) que comprovam presença (qualquer um basta)
    binaries: List[str]
    # teste de versão (lista de argv) — usado para confirmar que roda
    version_cmd: Optional[List[str]] = None
    # instaladores em ordem de prioridade; cada um retorna True/False
    installers: List[Callable[["InstallCtx"], bool]] = field(default_factory=list)
    # bancos de dados associados (função de download), opcional
    db_downloader: Optional[Callable[["InstallCtx"], bool]] = None
    note: str = ""

    def present(self) -> bool:
        if not any(have(b) for b in self.binaries):
            return False
        if self.version_cmd:
            cp = run(self.version_cmd, quiet=True)
            return cp.returncode == 0 or bool(cp.stdout or cp.stderr)
        return True


@dataclass
class InstallCtx:
    db_dir: str
    download_dbs: bool


# ----------------------------------------------------------------------
# Helpers de instalação (os "aliases"/fontes alternativas)
# ----------------------------------------------------------------------
def conda_install(pkgs: List[str], channels: List[str] = None) -> Callable:
    channels = channels or ["conda-forge", "bioconda"]

    def _do(ctx: InstallCtx) -> bool:
        fe = conda_frontend()
        if not fe:
            print(f"  {WARN} mamba/conda não encontrado — pulando fonte conda")
            return False
        chan = []
        for c in channels:
            chan += ["-c", c]
        cmd = [fe, "install", "-y", *chan, *pkgs]
        cp = run(cmd)
        return cp.returncode == 0
    return _do


def pip_install(pkgs: List[str]) -> Callable:
    def _do(ctx: InstallCtx) -> bool:
        cmd = [sys.executable, "-m", "pip", "install", "--upgrade", *pkgs]
        cp = run(cmd)
        return cp.returncode == 0
    return _do


def git_install(repo: str, postinstall: List[List[str]] = None) -> Callable:
    """Clona o repo em ~/tools/<nome> e roda pip install -e . (fallback final)."""
    postinstall = postinstall or [[sys.executable, "-m", "pip", "install", "-e", "."]]

    def _do(ctx: InstallCtx) -> bool:
        if not have("git"):
            print(f"  {WARN} git não encontrado — pulando fonte git")
            return False
        base = os.path.expanduser("~/tools")
        os.makedirs(base, exist_ok=True)
        name = repo.rstrip("/").split("/")[-1].replace(".git", "")
        dest = os.path.join(base, name)
        if not os.path.isdir(dest):
            if run(["git", "clone", "--depth", "1", repo, dest]).returncode != 0:
                return False
        cwd = os.getcwd()
        try:
            os.chdir(dest)
            for step in postinstall:
                if run(step).returncode != 0:
                    return False
        finally:
            os.chdir(cwd)
        return True
    return _do


# ----------------------------------------------------------------------
# Downloaders de bancos de dados
# ----------------------------------------------------------------------
def db_genomad(ctx: InstallCtx) -> bool:
    target = os.path.join(ctx.db_dir, "genomad_db")
    if os.path.isdir(target):
        print(f"  {OK} geNomad DB já existe em {target}")
        return True
    return run(["genomad", "download-database", ctx.db_dir]).returncode == 0


def db_checkv(ctx: InstallCtx) -> bool:
    # checkv cria checkv-db-vX.Y dentro do db_dir
    existing = [d for d in os.listdir(ctx.db_dir)
                if d.startswith("checkv-db")] if os.path.isdir(ctx.db_dir) else []
    if existing:
        print(f"  {OK} CheckV DB já existe ({existing[0]})")
        return True
    return run(["checkv", "download_database", ctx.db_dir]).returncode == 0


def db_phabox(ctx: InstallCtx) -> bool:
    target = os.path.join(ctx.db_dir, "phabox_db_v2")
    if os.path.isdir(target):
        print(f"  {OK} PhaBox2 DB já existe em {target}")
        return True
    # PhaBox2 distribui o DB via Zenodo (link muda por versão) -> ver --check
    print(f"  {WARN} Baixe o DB do PhaBox2 manualmente (Zenodo) e descompacte em {target}")
    print(f"        https://github.com/KennthShang/PhaBOX  (seção 'Download the database')")
    return os.path.isdir(target)


def db_vitap(ctx: InstallCtx) -> bool:
    target = os.path.join(ctx.db_dir, "ViTAP_db")
    if os.path.isdir(target):
        print(f"  {OK} ViTAP DB já existe em {target}")
        return True
    print(f"  {WARN} Baixe/gere o DB do ViTAP conforme o repo e aponte para {target}")
    print(f"        https://github.com/DamieFu/ViTAP")
    return os.path.isdir(target)


# ----------------------------------------------------------------------
# Catálogo de ferramentas
# ----------------------------------------------------------------------
def build_catalog() -> List[Tool]:
    return [
        Tool("geNomad", ["genomad"], ["genomad", "--version"],
             installers=[conda_install(["genomad"]),
                         pip_install(["genomad"]),
                         git_install("https://github.com/apcamargo/genomad")],
             db_downloader=db_genomad,
             note="identificação de contigs virais"),

        Tool("CheckV", ["checkv"], ["checkv", "--help"],
             installers=[conda_install(["checkv"]),
                         pip_install(["checkv"]),
                         git_install("https://bitbucket.org/berkeleylab/checkv")],
             db_downloader=db_checkv,
             note="QC de qualidade/genes virais"),

        Tool("vRhyme", ["vRhyme"], ["vRhyme", "--help"],
             installers=[conda_install(["vrhyme"]),
                         pip_install(["vRhyme"]),
                         git_install("https://github.com/AnantharamanLab/vRhyme")],
             note="binning viral"),

        Tool("VAMB", ["vamb"], ["vamb", "--help"],
             installers=[pip_install(["vamb"]),
                         conda_install(["vamb"]),
                         git_install("https://github.com/RasmussenLab/vamb")],
             note="binning (autoencoder)"),

        Tool("samtools", ["samtools"], ["samtools", "--version"],
             installers=[conda_install(["samtools"])],
             note="manipulação de BAM"),

        Tool("bwa-mem2", ["bwa-mem2"], ["bwa-mem2", "version"],
             installers=[conda_install(["bwa-mem2"])],
             note="mapeamento Illumina"),

        Tool("bowtie2", ["bowtie2"], ["bowtie2", "--version"],
             installers=[conda_install(["bowtie2"])],
             note="mapeamento Illumina (alternativa)"),

        Tool("minimap2", ["minimap2"], ["minimap2", "--version"],
             installers=[conda_install(["minimap2"])],
             note="mapeamento Nanopore"),

        Tool("CoverM", ["coverm"], ["coverm", "--version"],
             installers=[conda_install(["coverm"]),
                         git_install("https://github.com/wwood/CoverM")],
             note="abundância de vOTUs (metagenoma)"),

        Tool("Prodigal", ["prodigal", "prodigal-gv"], ["prodigal", "-v"],
             installers=[conda_install(["prodigal", "prodigal-gv"])],
             note="chamada de genes"),

        Tool("PhaBox2", ["phabox2", "phabox"], ["phabox2", "--help"],
             installers=[conda_install(["phabox"]),
                         pip_install(["phabox"]),
                         git_install("https://github.com/KennthShang/PhaBOX")],
             db_downloader=db_phabox,
             note="lifestyle (PhaTYP: temperate/lytic) + tax"),

        Tool("ViTAP", ["ViTAP", "vitap", "ViTAP.py"], None,
             installers=[git_install("https://github.com/DamieFu/ViTAP")],
             db_downloader=db_vitap,
             note="taxonomia viral"),
    ]


# ----------------------------------------------------------------------
# Lógica principal
# ----------------------------------------------------------------------
def preflight() -> int:
    """Verifica só o necessário para CONSTRUIR os envs isolados:
    conda/mamba, snakemake e git. Não checa as tools de bioinformática
    (elas vivem em envs isolados, não no PATH do base)."""
    print(f"{BOLD}== Preflight (orquestração) =={RESET}")
    fe = conda_frontend()
    snak = have("snakemake")
    git = have("git")
    print(f"  {(OK if fe else FAIL)} conda/mamba   {('('+fe+')') if fe else '— instale miniforge/miniconda'}")
    print(f"  {(OK if snak else FAIL)} snakemake     {'' if snak else '— mamba install -n base -c bioconda snakemake-minimal'}")
    print(f"  {(OK if git else WARN)} git           {'' if git else '— recomendado p/ fontes de fallback'}")
    if fe and snak:
        print(f"\n{GREEN}Preflight OK — pronto para construir os envs.{RESET}")
        return 0
    print(f"\n{RED}Preflight falhou.{RESET} Resolva os itens acima e tente de novo.")
    return 1


def try_update(tool: Tool) -> None:
    """Tenta atualizar uma ferramenta já instalada (conda/pip/git pull)."""
    fe = conda_frontend()
    print(f"  ↻ atualizando {tool.name}...")
    # 1) conda update (se houver frontend)
    if fe:
        # tenta o nome em minúsculo como pacote conda
        run([fe, "update", "-y", "-c", "conda-forge", "-c", "bioconda",
             tool.name.lower()], quiet=True)
    # 2) pip upgrade (silencioso; só funciona se for pacote pip)
    run([sys.executable, "-m", "pip", "install", "--upgrade", tool.name.lower()],
        quiet=True)
    # 3) git pull se foi instalado via repo em ~/tools
    repo_dir = os.path.expanduser(f"~/tools/{tool.name}")
    if os.path.isdir(os.path.join(repo_dir, ".git")):
        cwd = os.getcwd()
        try:
            os.chdir(repo_dir)
            run(["git", "pull", "--ff-only"], quiet=True)
        finally:
            os.chdir(cwd)


def try_install(tool: Tool, ctx: InstallCtx) -> bool:
    print(f"\n{BOLD}Instalando {tool.name}{RESET} ({tool.note})")
    for i, installer in enumerate(tool.installers, 1):
        print(f"  fonte {i}/{len(tool.installers)}...")
        if installer(ctx):
            if tool.present():
                print(f"  {OK} {tool.name} instalado com sucesso")
                return True
            print(f"  {WARN} comando rodou mas {tool.name} ainda não detectado")
        else:
            print(f"  {WARN} fonte {i} falhou — tentando próxima")
    print(f"  {FAIL} não consegui instalar {tool.name} automaticamente")
    return False


def main():
    ap = argparse.ArgumentParser(description="Checa/instala as tools do pipeline viral")
    ap.add_argument("--check", action="store_true", help="apenas relatório")
    ap.add_argument("--install", action="store_true", help="instala o que faltar")
    ap.add_argument("--update", action="store_true", help="atualiza as tools já instaladas")
    ap.add_argument("--download-dbs", action="store_true", help="baixa os bancos de dados")
    ap.add_argument("--db-dir", default=os.path.expanduser("~/viral_dbs"),
                    help="diretório dos bancos de dados")
    ap.add_argument("--only", default=None,
                    help="restringe a checagem/instalação a estas tools (lista por vírgula)")
    ap.add_argument("--preflight", action="store_true",
                    help="checa só conda/mamba + snakemake + git (p/ envs isolados)")
    args = ap.parse_args()

    if args.preflight:
        sys.exit(preflight())

    if not (args.check or args.install):
        args.check = True

    os.makedirs(args.db_dir, exist_ok=True)
    ctx = InstallCtx(db_dir=os.path.abspath(args.db_dir),
                     download_dbs=args.download_dbs)

    fe = conda_frontend()
    print(f"{BOLD}== Pipeline viral — checagem de ferramentas =={RESET}")
    print(f"conda/mamba: {OK + ' ' + fe if fe else FAIL + ' nenhum (instale miniforge)'}")
    print(f"db-dir     : {ctx.db_dir}\n")

    catalog = build_catalog()
    if args.only:
        wanted = {n.strip().lower() for n in args.only.split(",")}
        catalog = [t for t in catalog if t.name.lower() in wanted]
        if not catalog:
            print(f"{FAIL} --only não casou com nenhuma tool conhecida")
            sys.exit(2)

    status = {}     # nome -> (presente?, instalado_agora?)

    for tool in catalog:
        present = tool.present()
        installed_now = False
        if present:
            print(f"{OK} {tool.name:10s} — presente   ({tool.note})")
            if args.update:
                try_update(tool)
        else:
            print(f"{FAIL} {tool.name:10s} — FALTANDO   ({tool.note})")
            if args.install:
                installed_now = try_install(tool, ctx)
                present = installed_now
        status[tool.name] = present

        # bancos de dados
        if present and tool.db_downloader and args.download_dbs:
            print(f"  ↳ DB de {tool.name}:")
            tool.db_downloader(ctx)

    # ------------------------------------------------------------------
    # Relatório final
    # ------------------------------------------------------------------
    print(f"\n{BOLD}================= RESUMO ================={RESET}")
    missing = [n for n, ok in status.items() if not ok]
    for tool in catalog:
        mark = OK if status[tool.name] else FAIL
        print(f"  {mark} {tool.name}")

    if missing:
        print(f"\n{RED}{BOLD}FALTANDO:{RESET} {', '.join(missing)}")
        print("Sugestões:")
        print("  • rode novamente com --install")
        print("  • PhaBox2/ViTAP às vezes precisam de DB manual (ver mensagens acima)")
        print("  • confira se mamba/conda está ativo: `conda activate base`")
        sys.exit(1)
    else:
        print(f"\n{GREEN}{BOLD}Tudo pronto!{RESET} Todas as ferramentas estão disponíveis.")
        if not args.download_dbs:
            print("Lembre de baixar os bancos: --download-dbs --db-dir <dir>")
        sys.exit(0)


if __name__ == "__main__":
    main()
