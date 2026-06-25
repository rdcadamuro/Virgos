#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phages.py — launcher amigável do pipeline viral.

Uso interativo (recomendado):
    python phages.py
  -> escolhe o modo (genome/metagenome), aponta a entrada e a saída,
     e seleciona o bloco pelo NÚMERO.

Uso direto (automação):
    # genoma: pasta cheia de genomas completos
    python phages.py --mode genome --input-folder /dados/genomas -o /saida --block 1

    # metagenoma: tabela de amostras (assembly + reads)
    python phages.py --mode metagenome --samples meus_samples.tsv -o /saida --block 1

Blocos:
    0  Checar / instalar / ATUALIZAR ferramentas + bancos de dados
    1  RODAR TUDO  (geNomad -> ... -> summary)         [recomendado]
    2  01 geNomad
    3  02 CheckV
    4  03 Binning + dereplicação (metagenoma) | filtro de contigs (genoma)
    5  04 PhaBox2  (lifestyle: temperate/lytic)
    6  05 ViTAP    (taxonomia)
    7  06 Abundância de vOTUs                            [só metagenoma]
    8  Summary final
"""
import argparse
import glob
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
SNAKEFILE = os.path.join(ROOT, "workflow", "Snakefile")
CHECK_TOOLS = os.path.join(ROOT, "workflow", "scripts", "check_tools.py")
FASTA_EXT = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")

# bloco -> (rótulo, alvo snakemake, modos onde aparece)
BLOCKS = [
    (0, "Checar / instalar / ATUALIZAR ferramentas + bancos", "__tools__", ("genome", "metagenome")),
    (1, "RODAR TUDO (recomendado)",                            "stage_summary", ("genome", "metagenome")),
    (2, "01 geNomad",                                          "stage_genomad", ("genome", "metagenome")),
    (3, "02 CheckV",                                           "stage_checkv",  ("genome", "metagenome")),
    (4, "03 Binning + dereplicação / filtro de contigs",       "stage_bins",    ("genome", "metagenome")),
    (5, "04 PhaBox2 (lifestyle temperate/lytic)",             "stage_phabox",  ("genome", "metagenome")),
    (6, "05 ViTAP (taxonomia)",                                "stage_vitap",   ("genome", "metagenome")),
    (7, "06 Abundância de vOTUs",                              "stage_abundance", ("metagenome",)),
    (8, "Summary final",                                       "stage_summary", ("genome", "metagenome")),
]


def ask(prompt, default=None):
    sfx = f" [{default}]" if default else ""
    r = input(f"{prompt}{sfx}: ").strip()
    return r or (default or "")


def find_fastas(folder):
    out = []
    for f in sorted(os.listdir(folder)):
        if f.lower().endswith(FASTA_EXT):
            out.append(os.path.join(folder, f))
    return out


def sample_name(path):
    base = os.path.basename(path)
    for ext in sorted(FASTA_EXT, key=len, reverse=True):
        if base.lower().endswith(ext):
            return base[: -len(ext)]
    return os.path.splitext(base)[0]


def build_samples_tsv(args, out_path):
    """Gera o samples.tsv conforme o modo/entrada. Retorna nº de amostras."""
    rows = []
    if args.mode == "genome":
        # genome também pode binar com vRhyme+VAMB se houver reads:
        #   --samples (tabela c/ reads)  ou  --input + --reads1 [--reads2] --read-type
        if args.samples:
            return args.samples
        elif args.input_folder:
            files = find_fastas(args.input_folder)
            if not files:
                sys.exit(f"Nenhum FASTA encontrado em {args.input_folder}")
            for f in files:
                rows.append((sample_name(f), os.path.abspath(f), "", "", "na"))
        elif args.input:
            if args.reads1:   # genoma de isolado COM reads -> bina
                rt = args.read_type or "ilu-paired"
                rows.append((sample_name(args.input), os.path.abspath(args.input),
                             os.path.abspath(args.reads1),
                             os.path.abspath(args.reads2) if args.reads2 else "", rt))
            else:             # genoma sem reads -> usa contigs virais direto
                rows.append((sample_name(args.input), os.path.abspath(args.input),
                             "", "", "na"))
        else:
            sys.exit("modo genome: use --input/--input-folder, ou --samples p/ binar com reads")

    else:  # metagenome
        if args.samples:
            # usa a tabela do usuário diretamente
            return args.samples
        elif args.input and args.reads1:
            rt = args.read_type or "ilu-paired"
            rows.append((sample_name(args.input), os.path.abspath(args.input),
                         os.path.abspath(args.reads1),
                         os.path.abspath(args.reads2) if args.reads2 else "", rt))
        else:
            sys.exit("modo metagenome: use --samples <tsv>  (ou --input + --reads1 "
                     "[--reads2] --read-type ilu-paired|ilu-single|nano)")

    with open(out_path, "w") as o:
        o.write("sample\tassembly\treads1\treads2\tread_type\n")
        for r in rows:
            o.write("\t".join(r) + "\n")
    return out_path


def run_block0(threads, db_dir=None):
    """Setup robusto: preflight -> envs isolados por ferramenta -> bancos."""
    db_dir = db_dir or ask("Diretório dos bancos de dados",
                            os.path.expanduser("~/viral_dbs"))

    # 1) preflight (conda/mamba + snakemake + git)
    print("\n--- 1/3 preflight ---")
    if subprocess.call([sys.executable, CHECK_TOOLS, "--preflight"]) != 0:
        return 1

    # 2) constrói UM env isolado por ferramenta (sem rodar nada)
    print("\n--- 2/3 construindo envs isolados (pode demorar) ---")
    rc = subprocess.call(
        ["snakemake", "-s", SNAKEFILE, "--use-conda", "--conda-create-envs-only",
         "--cores", str(threads), "setup_envs", "--config", f"db_dir={db_dir}"],
        cwd=ROOT)
    if rc != 0:
        print("Falha ao construir os envs."); return rc

    # 3) baixa os bancos (cada um no env da própria ferramenta)
    print("\n--- 3/3 baixando bancos de dados ---")
    rc = subprocess.call(
        ["snakemake", "-s", SNAKEFILE, "--use-conda", "--cores", str(threads),
         "setup_dbs", "--config", f"db_dir={db_dir}"],
        cwd=ROOT)
    if rc == 0:
        print(f"\n✔ Setup concluído. Bancos em: {db_dir}")
        print("  Agora ajuste os caminhos em config/config.yaml -> databases:")
    return rc


def run_snakemake(target, mode, outdir, samples_path, threads, dry, extra):
    cmd = ["snakemake", "-s", SNAKEFILE, "--use-conda",
           "--cores", str(threads), target,
           "--config", f"mode={mode}", f"outdir={os.path.abspath(outdir)}",
           f"samples={os.path.abspath(samples_path)}"]
    if dry:
        cmd.append("-n")
    cmd += extra
    print("\n>>", " ".join(cmd), "\n")
    return subprocess.call(cmd, cwd=ROOT)


def choose_interactive(args):
    print("=" * 60)
    print("  phages.py — pipeline viral (geNomad -> ... -> ecologia)")
    print("=" * 60)
    if not args.mode:
        m = ask("Modo  (1=genome  2=metagenome)", "2")
        args.mode = "genome" if m.strip() in ("1", "genome", "g") else "metagenome"

    if not (args.input or args.input_folder or args.samples):
        if args.mode == "genome":
            p = ask("Pasta com genomas completos (ou 1 arquivo FASTA)")
            if os.path.isdir(p):
                args.input_folder = p
            else:
                args.input = p
        else:
            args.samples = ask("Tabela de amostras (samples.tsv: assembly+reads+read_type)")

    if not args.output:
        args.output = ask("Pasta de saída (-o)", os.path.abspath("results"))
    if not args.threads:
        args.threads = int(ask("Núcleos (threads)", "8"))

    print("\nBlocos disponíveis:")
    for num, label, _t, modes in BLOCKS:
        if args.mode in modes:
            print(f"   {num}) {label}")
    sel = ask("\nDigite o número do bloco", "1")
    return int(sel)


def main():
    ap = argparse.ArgumentParser(description="Launcher do pipeline viral")
    ap.add_argument("--mode", choices=["genome", "metagenome"])
    ap.add_argument("--input", help="1 arquivo FASTA (genoma único / assembly único)")
    ap.add_argument("--input-folder", help="pasta cheia de FASTA (modo genome)")
    ap.add_argument("--samples", help="samples.tsv pronto (modo metagenome)")
    ap.add_argument("--reads1"); ap.add_argument("--reads2")
    ap.add_argument("--read-type", choices=["ilu-paired", "ilu-single", "nano"])
    ap.add_argument("-o", "--output", help="pasta de saída")
    ap.add_argument("-t", "--threads", type=int)
    ap.add_argument("--db-dir", help="diretório dos bancos de dados (bloco 0)")
    ap.add_argument("--block", type=int, help="número do bloco (0-8); omitido = menu")
    ap.add_argument("-n", "--dry-run", action="store_true", help="só mostra o DAG")
    ap.add_argument("extra", nargs="*", help="args extras repassados ao snakemake")
    args = ap.parse_args()

    block = args.block
    if block is None:
        block = choose_interactive(args)

    # defaults se veio tudo por flag
    args.output = args.output or os.path.abspath("results")
    args.threads = args.threads or 8

    # bloco 0 não precisa de amostras
    if block == 0:
        sys.exit(run_block0(args.threads, args.db_dir))

    if not args.mode:
        sys.exit("Defina --mode genome|metagenome")

    blk = next((b for b in BLOCKS if b[0] == block), None)
    if not blk:
        sys.exit(f"Bloco inválido: {block}")
    if args.mode not in blk[3]:
        sys.exit(f"Bloco {block} ({blk[1]}) não se aplica ao modo {args.mode}")

    os.makedirs(args.output, exist_ok=True)
    samples_path = build_samples_tsv(args, os.path.join(args.output, "samples.tsv"))

    rc = run_snakemake(blk[2], args.mode, args.output, samples_path,
                       args.threads, args.dry_run, args.extra)
    if rc == 0:
        print(f"\n✔ Bloco {block} concluído. Saída em: {os.path.abspath(args.output)}")
    sys.exit(rc)


if __name__ == "__main__":
    main()
