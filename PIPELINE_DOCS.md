# Documentação Técnica — phages viral pipeline

> **Status:** validado end-to-end em 2026-06-24 com amostra de genoma completo (MS13934, Caudoviricetes).
> Inclui: correção da cadeia de falha silenciosa do vRhyme (§10.8), CheckV no bin inteiro / N-linking (§12)
> e verificação inteligente anti-erro-silencioso (§11).

---

## Índice

1. [Visão geral](#1-visão-geral)
2. [Arquitetura e fluxo de trabalho](#2-arquitetura-e-fluxo-de-trabalho)
3. [Ferramentas e versões](#3-ferramentas-e-versões)
4. [Bancos de dados](#4-bancos-de-dados)
5. [Requisitos do sistema](#5-requisitos-do-sistema)
6. [Instalação passo a passo (WSL2)](#6-instalação-passo-a-passo-wsl2)
7. [Configuração](#7-configuração)
8. [Uso](#8-uso)
9. [Descrição das saídas](#9-descrição-das-saídas)
10. [Problemas conhecidos e correções aplicadas](#10-problemas-conhecidos-e-correções-aplicadas)
11. [Verificação de integridade da execução (anti-erro-silencioso)](#11-verificação-de-integridade-da-execução-anti-erro-silencioso)
12. [CheckV no bin inteiro — N-linking (metodologia vRhyme)](#12-checkv-no-bin-inteiro--n-linking-metodologia-vrhyme)
13. [Anotação funcional: resistoma + mobiloma + enriquecimento](#13-anotação-funcional-resistoma--mobiloma--enriquecimento)
14. [Submissão ao ENA: pacote MIUViG dos uViGs](#14-submissão-ao-ena-pacote-miuvig-dos-uvigs)
15. [Resultado da validação — amostra MS13934](#15-resultado-da-validação--amostra-ms13934)
16. [Referências bibliográficas (ferramentas a citar)](#16-referências-bibliográficas-ferramentas-a-citar)

---

## 1. Visão geral

Pipeline Snakemake para identificação, refinamento, anotação e (no modo metagenoma) quantificação de fagos a partir de sequências metagenômicas ou de genomas completos.

**Dois modos de operação:**

| Modo | Entrada | Etapas ativas |
|------|---------|---------------|
| `genome` | FASTA de genoma/assembly (com ou sem reads) | geNomad → CheckV → [vRhyme+VAMB se houver reads] → derep → PhaBox2 → VITAP → summary |
| `metagenome` | assembly + reads (Illumina ou Nanopore) | geNomad → CheckV → vRhyme+VAMB → derep → PhaBox2 → VITAP → **abundância** → summary |

O binning (vRhyme+VAMB) é disparado pela **presença de reads**, não pelo modo. O modo `metagenome` adiciona exclusivamente a etapa de abundância de vOTUs.

---

## 2. Arquitetura e fluxo de trabalho

```
ENTRADA: FASTA [+ reads]
         │
         ▼
┌─────────────────────────────────────────────────────────────────────┐
│  BLOCO 01 — geNomad                                                 │
│  Identificação de sequências virais por modelos de linguagem + MMseqs2│
│  Saída: contigs virais (score ≥ 0.7) → 01_genomad/viral_contigs/   │
└───────────────────────────┬─────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────────┐
│  BLOCO 02 — CheckV                                                  │
│  QC: genes virais, completude, contaminação                         │
│  Saída: quality_summary.tsv → 02_checkv/                           │
└───────────────────────────┬─────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────────┐
│  BLOCO 03 — Binning + Dereplicação                                  │
│  ┌─────────────────────┐  ┌──────────────────────────┐             │
│  │  vRhyme             │  │  VAMB                    │             │
│  │  (co-variância reads│  │  (variational autoencoder│             │
│  │   → bins de fagos)  │  │   → bins alternativos)   │             │
│  └──────────┬──────────┘  └─────────────┬────────────┘             │
│             └──────────────┬────────────┘                          │
│                            ▼                                        │
│     link_bins: scaffolds do bin unidos por N's (§12)               │
│                            ▼                                        │
│     checkv_bins: CheckV no BIN inteiro (completude do vMAG)         │
│                            ▼                                        │
│             filter_units (min 1 gene viral CheckV)                  │
│                            ▼                                        │
│             derep ANI ≥ 95% / AF ≥ 85% → votus.fna                │
│  Saída: 03_bins_vrhyme_vamb/  +  02_checkv/<amostra>_bins/         │
└───────────────────────────┬─────────────────────────────────────────┘
                            │
               ┌────────────┴────────────┐
               ▼                         ▼ (só metagenoma)
┌──────────────────────────┐  ┌─────────────────────────────────────┐
│  BLOCO 04 — PhaBox2      │  │  BLOCO 06 — Abundância (CoverM)     │
│  PhaMer: identificação   │  │  Mapeia reads → vOTUs (bwa-mem2     │
│  PhaTYP: lifestyle       │  │  / bowtie2 / minimap2)              │
│  PhaGCN: genus cluster   │  │  Saída: votu_abundance_tpm.tsv      │
│  CHERRY: host prediction │  │  (TPM, meandepth, covfraction)      │
│  Saída: 04_phabox2/      │  └──────────────┬──────────────────────┘
└──────────────┬───────────┘                 │
               ▼                             │
┌──────────────────────────┐                 │
│  BLOCO 05 — VITAP v1.10  │                 │
│  Taxonomia: BLAST vs ICTV│                 │
│  + grafo bipartido multi │                 │
│  -taxon (Reino→Espécie)  │                 │
│  Saída: 05_vitap/        │                 │
└──────────────┬───────────┘                 │
               └──────────────┬──────────────┘
                              ▼
┌─────────────────────────────────────────────────────────────────────┐
│  BLOCO FINAL — Tabela integrada                                     │
│  viral_bins_table.tsv + per_sample_summary.tsv                     │
│  Saída: 06_summary_final/ (genome) ou 07_summary_final/ (metag.)   │
└───────────────────────────┬─────────────────────────────────────────┘
                            ▼
┌─────────────────────────────────────────────────────────────────────┐
│  healthcheck — verificação inteligente anti-erro-silencioso (§11)   │
│  varre logs + cruza consistência entre etapas → RUN_HEALTHCHECK.txt │
└─────────────────────────────────────────────────────────────────────┘
```

O Snakemake gerencia o DAG de dependências: cada regra só dispara quando seus
inputs estão prontos, e resultados intermediários são reutilizados em re-runs.
Cada ferramenta roda em seu **próprio ambiente conda isolado** (`--use-conda`),
eliminando conflitos de dependências.

---

## 3. Ferramentas e versões

| Bloco | Ferramenta | Versão mínima | Conda env | Referência |
|-------|-----------|--------------|-----------|-----------|
| 01 | **geNomad** | ≥ 1.7 | `genomad` | Camargo et al., 2023 |
| 02 | **CheckV** | ≥ 1.0.1 | `checkv` | Nayfach et al., 2021 |
| 02 | **Diamond** | ≥ 2.1 | `checkv` | Buchfink et al., 2021 |
| 02 | **HMMER** | ≥ 3.3 | `checkv` | Eddy, 2011 |
| 02 | **Prodigal-GV** | — | `checkv` | Camargo et al., 2023 |
| 03 | **vRhyme** | v1.1.0 (env c/ `python<3.12`, `scikit-learn<1.3`, `setuptools<81` — ver §10.8) | `vrhyme` | Kieft et al., 2022 |
| 03 | **VAMB** | ≥ 4.0 | `vamb` | Nissen et al., 2021 |
| 03 | **BLAST+** | ≥ 2.13 | `derep` | Camacho et al., 2009 |
| 03 | **bwa-mem2** | — | `mapping` | Vasimuddin et al., 2019 |
| 03 | **bowtie2** | — | `mapping` | Langmead & Salzberg, 2012 |
| 03 | **minimap2** | — | `mapping` | Li, 2018 |
| 03 | **samtools** | ≥ 1.17 | `mapping` | Li et al., 2009 |
| 04 | **PhaBox2** | ≥ 2.0 | `phabox` | Shang et al., 2024 |
| 05 | **VITAP** | = 1.10 | `vitap` | Zheng et al., 2025 |
| 06 | **CoverM** | ≥ 0.7 | `abundance` | — |
| orq.| **Snakemake** | ≥ 9.0 | base/viralpipe | Mölder et al., 2021 |

> **Notas de compatibilidade (críticas):**
> - **vRhyme** falha em silêncio sem os 3 pins do env (`python<3.12`,
>   `scikit-learn<1.3`, `setuptools<81`) — ver [§10.8](#108-vrhyme-falha-silenciosa-em-cadeia-crítico).
> - **PhaBox2** v2.1.13 requer patch para pandas ≥ 2.0, re-aplicado
>   automaticamente — ver [§10.3](#103-phabox2-v2113-incompatibilidade-com-pandas--20).

---

## 4. Bancos de dados

| Ferramenta | Banco | Tamanho aprox. | Localização padrão |
|-----------|-------|---------------|--------------------|
| geNomad | `genomad_db` | ~3 GB | `~/viral_dbs/genomad_db/` |
| CheckV | `checkv-db-v1.5` | ~2 GB | `~/viral_dbs/checkv-db-v1.5/` |
| PhaBox2 | `phabox_db_v2` | ~8 GB | `~/viral_dbs/phabox_db_v2/` |
| VITAP | `DB_hybrid_MSL37_RefSeq209_IMGVR` | ~2.3 GB | `~/viral_dbs/ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR/` |

**VITAP — atenção importante:**
- O banco é distribuído no Figshare como `.zip`. Após extração, cria uma subpasta
  `DB_hybrid_MSL37_RefSeq209_IMGVR/` dentro do destino — o `config.yaml` deve
  apontar para **essa subpasta**, não para o diretório pai.
- O arquivo `VMR_genome_hybrid_RefSeqs_IMGVR.dmnd` distribuído no Figshare foi
  compilado com uma versão antiga do Diamond. É necessário **reconstruir** com
  a versão instalada:
  ```bash
  conda activate <env_vitap>
  diamond makedb \
    --in ~/viral_dbs/ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR/VMR_genome_hybrid_RefSeqs_IMGVR.faa \
    --db ~/viral_dbs/ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR/VMR_genome_hybrid_RefSeqs_IMGVR \
    --threads 4
  ```
- O `uniref90.dmnd` (banco UniRef90 completo, ~20 GB) **não** está incluso.
  O pipeline cria automaticamente um arquivo de fallback vazio
  (`target_uniref90_taxa_fallback.out`) quando o arquivo não existe, pulando
  essa etapa de refinamento de linhagem sem comprometer a saída principal.

---

## 5. Requisitos do sistema

| Componente | Requisito |
|------------|----------|
| SO | Linux (testado: Ubuntu 22.04 / WSL2) |
| RAM | ≥ 8 GB (16 GB recomendado para geNomad sem `--splits`) |
| Disco | ≥ 20 GB para bancos de dados + espaço para resultados |
| CPU | ≥ 4 cores (8 recomendado) |
| Python | ≥ 3.10 (orquestrador); envs isolados por ferramenta |
| Conda | miniconda3 ou miniforge3 |

**WSL2 — restrições específicas:**
- O diretório de saída **deve estar no filesystem Linux** (`/home/...`), não em
  `/mnt/c/...` (NTFS). O MMseqs2 (usado pelo geNomad) falha ao alocar memória
  mapeada em NTFS.
- Com 8 GB de RAM, o geNomad requer `--splits 4` (ou superior) em seu parâmetro
  `extra` no `config.yaml` para evitar OOM no índice MMseqs2.
- Ao invocar o pipeline via `wsl bash -c "..."`, o conda deve ser
  explicitamente inicializado:
  ```bash
  source /home/rafa/miniconda3/etc/profile.d/conda.sh
  conda activate base
  python3 phages.py ...
  ```

---

## 6. Instalação passo a passo (WSL2)

### 6.1 Conda / Miniforge

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
~/miniforge3/bin/conda init bash
exec bash
```

### 6.2 Ambiente de orquestração

```bash
mamba create -n viralpipe -c conda-forge -c bioconda snakemake-minimal pandas
conda activate viralpipe
```

> O `snakemake-minimal` instala o Snakemake sem dependências pesadas de UI.
> O `pandas` é necessário para o Snakefile ler o `samples.tsv`.

### 6.3 Clonar o repositório

```bash
git clone https://github.com/<usuario>/phages-pipeline.git ~/viral_pipeline
cd ~/viral_pipeline
```

### 6.4 Baixar os bancos de dados (Bloco 0)

```bash
conda activate viralpipe
python phages.py --block 0
```

O Bloco 0 cria todos os ambientes conda por ferramenta e baixa os bancos
automaticamente onde possível. Para PhaBox2 e VITAP, exibe as URLs de download
manual quando necessário.

**Após baixar o banco VITAP**, reconstruir o arquivo Diamond (ver seção 4).

### 6.5 Editar o config.yaml

Atualizar os caminhos em `config/config.yaml` conforme a localização dos bancos:

```yaml
databases:
  genomad: "/home/<user>/viral_dbs/genomad_db"
  checkv:  "/home/<user>/viral_dbs/checkv-db-v1.5"
  phabox:  "/home/<user>/viral_dbs/phabox_db_v2"
  vitap:   "/home/<user>/viral_dbs/ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR"

genomad:
  min_virus_score: 0.7
  extra: "--cleanup --splits 4"   # --splits 4 obrigatório em WSL2 com 8 GB RAM
```

---

## 7. Configuração

Arquivo principal: [`config/config.yaml`](config/config.yaml)

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `genomad.min_virus_score` | `0.7` | Score mínimo geNomad para reter contig |
| `genomad.extra` | `"--cleanup --splits 4"` | Flags extras para geNomad (obrigatório `--splits 4` em WSL2/8GB) |
| `mapping.illumina_mapper` | `"bwa-mem2"` | Mapeador para reads Illumina (`bwa-mem2` ou `bowtie2`) |
| `binning.min_contig_len` | `2000` | Tamanho mínimo de contig para binning |
| `binning.run_vrhyme` | `true` | Ativa vRhyme |
| `binning.run_vamb` | `true` | Ativa VAMB |
| `checkv.min_viral_genes` | `1` | Mínimo de genes virais CheckV para reter bin |
| `derep.min_ani` | `95` | ANI mínima para dereplicação de vOTUs (%) |
| `derep.min_af` | `85` | Alignment fraction mínima para dereplicação (%) |
| `phabox.task` | `"end_to_end"` | Modo PhaBox2 (roda PhaMer+PhaGCN+PhaTYP+CHERRY) |
| `abundance.min_read_pct_id` | `95` | Identidade mínima do read para contagem de abundância |
| `abundance.min_read_aln_pct` | `75` | Fração mínima do read alinhada |
| `threads.*` | `8` | Threads por etapa (escalado para `--cores` do Snakemake) |

---

## 8. Uso

### Modo interativo

```bash
conda activate viralpipe
cd ~/viral_pipeline
python phages.py
```

### Modo direto — exemplos

```bash
# Um genoma completo com reads paired-end
python phages.py \
  --mode genome \
  --input /dados/genoma.fasta \
  --reads1 /dados/R1.fastq.gz \
  --reads2 /dados/R2.fastq.gz \
  --read-type ilu-paired \
  -o /resultados/amostra1 \
  --block 1

# Pasta com múltiplos genomas (sem reads)
python phages.py \
  --mode genome \
  --input-folder /dados/genomas/ \
  -o /resultados/lote1 \
  --block 1

# Metagenoma — uma amostra
python phages.py \
  --mode metagenome \
  --input /dados/assembly.fasta \
  --reads1 /dados/R1.fastq.gz \
  --reads2 /dados/R2.fastq.gz \
  --read-type ilu-paired \
  -o /resultados/meta1 \
  --block 1

# Metagenoma — múltiplas amostras via tabela
python phages.py \
  --mode metagenome \
  --samples config/samples.tsv \
  -o /resultados/multi \
  --block 1

# Dry-run (mostra DAG sem executar)
python phages.py --mode genome --input genoma.fasta -o /tmp/test --block 1 -- -n

# Executar apenas uma etapa específica (ex: só PhaBox2)
python phages.py --mode genome --input genoma.fasta -o /resultados -o /res --block 5
```

### Tabela de amostras (`samples.tsv`)

Para múltiplas amostras no modo metagenoma:

```tsv
sample	assembly	reads1	reads2	read_type
AMOSTRA_A	/dados/A.fasta	/dados/A_R1.fq.gz	/dados/A_R2.fq.gz	ilu-paired
AMOSTRA_B	/dados/B.fasta	/dados/B_R1.fq.gz		ilu-single
AMOSTRA_C	/dados/C.fasta	/dados/C_nano.fq.gz		nano
```

Tipos de read válidos: `ilu-paired`, `ilu-single`, `nano`.

---

## 9. Descrição das saídas

### Estrutura de diretórios

```
<saida>/
├── 01_genomad/
│   ├── <amostra>/                  logs geNomad por etapa
│   └── viral_contigs/
│       ├── <amostra>.viral.fna     contigs virais filtrados
│       └── <amostra>.keep_ids.txt  IDs retidos (score ≥ min_virus_score)
│
├── 02_checkv/
│   ├── <amostra>/                  [sem reads] CheckV por contig
│   │   ├── quality_summary.tsv     completude, genes virais, qualidade
│   │   ├── completeness.tsv
│   │   ├── contamination.tsv
│   │   ├── viruses.fna
│   │   └── proviruses.fna
│   └── <amostra>_bins/             [com reads] CheckV no BIN inteiro (N-linking, §12)
│       └── quality_summary.tsv     completude/qualidade por vMAG (bin)
│
├── 03_bins_vrhyme_vamb/
│   ├── mapping/
│   │   └── <amostra>.sorted.bam    reads mapeados (para binning por co-cobertura)
│   ├── <amostra>/
│   │   ├── linked_bins.fna         bins com scaffolds unidos por N's (p/ CheckV, §12)
│   │   ├── passed.fna              bins que passaram no filtro de genes virais
│   │   └── passed_units.tsv        tabela de bins aprovados
│   ├── collected/
│   │   └── <amostra>.bins.tsv      metadados dos bins coletados
│   └── derep/
│       ├── votus.fna               representantes dereplicados (vOTUs)
│       ├── clusters.tsv            clusters de dereplicação
│       ├── contig_to_votu.tsv      mapa contig → vOTU representante
│       └── all_passed_units.tsv    todos os bins de todas as amostras
│
├── 04_phabox2/
│   └── final_prediction/
│       ├── final_prediction_summary.tsv   tabela completa integrada
│       ├── phamer_prediction.tsv          identificação como vírus
│       ├── phatyp_prediction.tsv          lifestyle (temperate/lytic)
│       ├── phagcn_prediction.tsv          cluster de gênero
│       ├── cherry_prediction.tsv          predição de hospedeiro
│       └── phavip_prediction.tsv          virulência
│
├── 05_vitap/
│   ├── vitap_taxonomy.tsv              taxonomia principal (cópia normalizada)
│   ├── best_determined_lineages.tsv    linhagens com maior confiança
│   ├── all_lineages.tsv                todas as linhagens candidatas
│   └── *_bipartite.graph               grafos por nível taxonômico
│
├── 06_votu_abundance/   [apenas metagenoma]
│   ├── votu_abundance_tpm.tsv          matriz vOTU × amostra em TPM
│   ├── votu_abundance_meandepth.tsv    profundidade média
│   └── votu_abundance_coveredfraction.tsv  fração coberta
│
├── 06_summary_final/    [modo genome]   ou   07_summary_final/  [modo metagenoma]
│   ├── viral_bins_table.tsv    ← RESULTADO PRINCIPAL
│   ├── per_sample_summary.tsv  ← RESUMO POR AMOSTRA
│   └── RUN_HEALTHCHECK.txt     ← VERIFICAÇÃO DE INTEGRIDADE (§11)
│
└── logs/                    logs individuais por regra e amostra
    ├── 01_genomad/  02_checkv/  02_checkv_bins/
    ├── 03_mapping/  03_vrhyme/  03_vamb/  03_derep.log
    └── 04_phabox.log  05_vitap.log  ...
```

### Colunas da tabela principal (`viral_bins_table.tsv`)

| Coluna | Descrição |
|--------|-----------|
| `votu_rep` | ID do vOTU representante (contig original) |
| `rep_unit` | bin de origem do representante |
| `rep_sample` | amostra de origem |
| `all_samples` | todas as amostras onde o vOTU foi detectado |
| `binner` | ferramenta de binning (`vamb` ou `vrhyme`) |
| `checkv_quality` | qualidade CheckV: Complete / High-quality / Medium-quality / Low-quality |
| `completeness` | completude estimada pelo CheckV (%) |
| `lifestyle` | `lysogenic`, `lytic` ou `-` (PhaTYP) |
| `lifestyle_score` | score de confiança do lifestyle |
| `taxonomy` | linhagem taxonômica VITAP (vazia se abaixo do threshold) |
| `n_members` | nº de bins no cluster do vOTU |
| `mean_tpm` | [metagenoma] abundância média em TPM |
| `prevalence` | [metagenoma] nº de amostras onde detectado |

### Colunas do resumo (`per_sample_summary.tsv`)

| Coluna | Descrição |
|--------|-----------|
| `sample` | nome da amostra |
| `n_votus_origin` | nº de vOTUs oriundos desta amostra |
| `n_lysogenic` | nº de vOTUs lisogênicos |
| `n_lytic` | nº de vOTUs líticos |
| `n_unknown` | nº sem lifestyle definido |
| `q_Complete` / `q_High-quality` / `q_Medium-quality` / `q_Low-quality` | nº de vOTUs por categoria de qualidade CheckV |

---

## 10. Problemas conhecidos e correções aplicadas

### 10.1 geNomad: OOM em WSL2 com 8 GB RAM

**Sintoma:** `Can not allocate entries memory in IndexTable::initMemory` — o MMseqs2
tenta indexar ~228 mil perfis de uma vez.

**Correção:** adicionar `--splits 4` ao parâmetro `genomad.extra` no `config.yaml`:
```yaml
genomad:
  extra: "--cleanup --splits 4"
```
Isso divide o índice em 4 partes menores, reduzindo o pico de memória.

### 10.2 geNomad / MMseqs2: falha em filesystem NTFS (`/mnt/c/`)

**Sintoma:** `IndexTable::initMemory` falha ao alocar memória mapeada em arquivo.

**Correção:** diretório de saída **deve estar em ext4** (ex: `/home/rafa/results/`).
Acessar de Windows via `\\wsl.localhost\Ubuntu-22.04\home\rafa\results\`.

### 10.3 PhaBox2 v2.1.13: incompatibilidade com pandas ≥ 2.0

**Sintoma:** `TypeError: Invalid value 'X' for dtype int64` ou `LossySetitemError`
em `phagcn.py`, `cherry.py`, `votu.py`, `phatyp.py`.

**Correção:** patch aplicado em três níveis:

**a) Shim global** em `phabox2.py` (após `import pandas as pd`):
```python
import pandas.core.internals.blocks as _pd_blk
_orig_coerce = _pd_blk.Block.coerce_to_target_dtype
def _patched_coerce(self, element, raise_on_upcast=False):
    return _orig_coerce(self, element, raise_on_upcast=False)
_pd_blk.Block.coerce_to_target_dtype = _patched_coerce
```

**b) Conversão explícita** no site específico de cada módulo:
- `phagcn.py`: `df["cluster"] = df["cluster"].astype(object)`
- `cherry.py`: `df["Score"] = df["Score"].astype(object)`
- `votu.py`: `cluster_df["cluster"] = cluster_df["cluster"].astype(object)`

**Durabilidade (importante):** esses patches vivem em
`.snakemake/conda/<hash>_/lib/python3.12/site-packages/phabox2/`, que é
**descartável** (gitignored, reconstruído a cada `--block 0`). Para não se perderem,
o script [`workflow/scripts/patch_phabox2.py`](workflow/scripts/patch_phabox2.py)
re-aplica os 4 patches de forma **idempotente** e é chamado automaticamente no
início da regra `phabox` (em `phabox.smk`). Assim a correção sobrevive a
reinstalações sem intervenção manual. (Pin de `pandas<2` não serve: o env do
PhaBox2 usa Python 3.12, e pandas 1.x não tem wheel para 3.12.)

### 10.4 VITAP: banco Diamond desatualizado (Figshare)

**Sintoma:** `Error: Database was built with an older version of Diamond` ao rodar
o BLAST contra o banco ICTV.

**Correção:** reconstruir o `.dmnd` com a versão instalada (Diamond 2.1.16):
```bash
diamond makedb \
  --in ~/viral_dbs/ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR/VMR_genome_hybrid_RefSeqs_IMGVR.faa \
  --db ~/viral_dbs/ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR/VMR_genome_hybrid_RefSeqs_IMGVR \
  --threads 4
# resultado: 725.659 sequências, 6,3 s
```

### 10.5 VITAP: `uniref90.dmnd` ausente (banco não incluído no Figshare)

**Sintoma:** `Error opening file .../uniref90.dmnd: No such file or directory`
na etapa "Filtering lineage information based on UniRef90".

**Contexto:** VITAP usa UniRef90 como fallback para confirmar a taxonomia quando
o BLAST ICTV não é conclusivo. O arquivo `.dmnd` (~20 GB) não é distribuído com
o banco principal.

**Correção implementada em `vitap.smk`:** antes de chamar o binário VITAP, criar
um arquivo de fallback vazio se o `.dmnd` não existir. O VITAP detecta o arquivo
existente e pula a etapa de BLAST UniRef90:
```bash
if [ ! -f {input.db}/uniref90.dmnd ]; then
    printf 'genome_id\ttaxa_name\tparticipation_index\ttaxon_level\n' \
        > {params.outdir}/target_uniref90_taxa_fallback.out
fi
```
**Impacto:** as linhagens são determinadas apenas pelo BLAST contra referências
ICTV. Para fragmentos com < 60% de completude, o resultado é frequentemente
sem atribuição taxonômica abaixo de Caudoviricetes, o que é biologicamente esperado.

### 10.6 VITAP: caminho do banco (`db_dir`) aponta para nível errado

**Sintoma:** `FileNotFoundError: No *.fasta found under db_dir`.

**Correção:** o Figshare ZIP extrai em `ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR/`.
O `config.yaml` deve apontar para a **subpasta**, não para `ViTAP_db/`:
```yaml
databases:
  vitap: "/home/rafa/viral_dbs/ViTAP_db/DB_hybrid_MSL37_RefSeq209_IMGVR"
```

### 10.7 WSL2: `snakemake` não encontrado ao invocar via `wsl bash -c`

**Sintoma:** `FileNotFoundError: [Errno 2] No such file or directory: 'snakemake'`

**Causa:** `wsl bash -c "..."` abre shell não-interativo sem ativar conda.

**Correção:** ativar conda explicitamente:
```bash
wsl -d "Ubuntu-22.04" bash -c "
source /home/rafa/miniconda3/etc/profile.d/conda.sh
conda activate base
python3 phages.py ...
"
```

### 10.8 vRhyme: falha silenciosa em cadeia (CRÍTICO)

**Sintoma (perigoso):** o vRhyme **não produzia bins**, mas o pipeline reportava
sucesso. A causa não era biológica — o vRhyme **quebrava no meio da execução** e a
regra mascarava o erro com `... || true`, então "0 bins" parecia legítimo quando na
verdade **o programa nunca terminava**. Os bins do teste vinham 100% do VAMB.

**Diagnóstico — três camadas encadeadas, todas escondidas pelo `|| true`:**

| Camada | Erro real | Causa |
|---|---|---|
| 1. Importação | `ModuleNotFoundError: No module named 'pkg_resources'` | `setuptools≥81` **removeu** `pkg_resources`; o vRhyme ainda o importa |
| 2. (falso positivo) | `mmseqs2 cannot be found` | artefato de testar `$ENV/bin/vRhyme` sem ativar o env (PATH). Na regra real o env é ativado — não é bug |
| 3. Classificação ML | `ValueError: node array from the pickle has an incompatible dtype` | os modelos ExtraTrees pré-treinados do vRhyme foram salvos com `scikit-learn<1.3`; o env tinha ≥1.3, que adicionou o campo `missing_go_to_left` na árvore |

**Correção (em `workflow/envs/vrhyme.yaml`):** três pins obrigatórios:
```yaml
dependencies:
  - "python<3.12"       # permite instalar scikit-learn<1.3 (sem wheel p/ 3.12+)
  - vrhyme
  - "scikit-learn<1.3"  # pickles do vRhyme usam o formato de árvore antigo
  - "setuptools<81"     # mantém pkg_resources disponível
  - samtools
```

**Correção estrutural (anti-mascaramento) em `binning.smk`:** o `|| true` foi
**removido** das regras `vrhyme` e `vamb`. Agora o código de saída é capturado e
classificado em três casos:
```bash
vRhyme ... ; rc=$?
if   [ "$rc" -eq 0 ];                              then echo "[vRhyme] OK ... Bins: N"
elif grep -qiE "no bins|not enough|no viral" LOG;  then echo "[vRhyme] OK -- 0 bins legítimo (sem co-cobertura)"
else echo "[vRhyme] *** ERRO REAL (exit=$rc) ***"; exit "$rc"   # ABORTA de propósito
fi
```
Assim "0 bins porque não havia sinal" (segue) fica **distinguível** de "0 bins
porque quebrou" (aborta). Mesma lógica aplicada ao VAMB.

**Lição:** padrões `comando || true` / `2>/dev/null` que engolem código de saída
são a principal fonte de falha silenciosa. A §11 descreve o verificador que caça
esse tipo de problema automaticamente.

---

## 11. Verificação de integridade da execução (anti-erro-silencioso)

Duas camadas garantem que **0 resultados ≠ erro escondido**:

### 11.1 Marcador de término por etapa (echo no log)

Cada regra escreve, no fim, uma linha explícita de status — distinguindo
"completou (exit 0)" de "0 resultados". Exemplos reais:
```
[geNomad] OK -- rodou ate o fim (exit 0). Sequencias virais: 5
[CheckV-bins] OK -- rodou ate o fim (exit 0). Bins (vMAGs) avaliados: 5
[vRhyme] OK (exit 0) -- rodou ate o fim. Bins formados: 0
[VITAP] OK -- rodou ate o fim (exit 0). Linhagens determinadas: 0
[VITAP] NOTA: 0 linhagens. VITAP COMPLETOU; ... NAO e erro.
```

### 11.2 Verificador inteligente consolidado

[`workflow/scripts/verify_run.py`](workflow/scripts/verify_run.py) roda
**automaticamente no fim** de toda execução (regra `healthcheck`, exigida por
`stage_summary`), gera `<saida>/0X_summary_final/RUN_HEALTHCHECK.txt` e imprime no
console. Detecta falha silenciosa por **dois mecanismos independentes**:

1. **Varredura de logs** — caça assinaturas fatais de alta precisão:
   `Traceback`, `No module named`, `command not found`, `CalledProcessError`,
   `Segmentation fault`, `core dumped`, `MemoryError`, `Cannot allocate`, `Killed`,
   e o marcador próprio `*** ERRO REAL`.
2. **Consistência entre etapas** — pega a falha que **não imprime erro nenhum**:
   ex. "geNomad achou N>0 mas CheckV avaliou 0 linhas" → `FALHA`; "N vOTUs mas 0
   predições do PhaBox2" → `FALHA`. (VITAP=0 linhagens é `AVISO`, não `FALHA`, pois
   é legítimo para fragmentos divergentes.)

Cada item recebe `OK` / `VAZIO-OK` (vazio legítimo) / `AVISO` / `FALHA`, com
veredito final: **INTEGRO** / **INTEGRO COM AVISOS** / **FALHA(S) DETECTADA(S)**.

Por padrão é **não-strict** (sempre completa e reporta alto, sem abortar). Para que
uma `FALHA` **aborte** o pipeline, rode o script com `--strict` (sai com código ≠0).

Uso avulso (a qualquer momento, sem re-rodar nada):
```bash
python workflow/scripts/verify_run.py \
    --outdir <saida> --mode genome --samples <saida>/samples.tsv \
    --out <saida>/06_summary_final/RUN_HEALTHCHECK.txt
```

**Exemplo de saída (run real MS13934):**
```
  [ OK ] 01 geNomad       5 contigs virais identificados
  [ OK ] 02 CheckV-bins   5 unidades avaliadas (completude/qualidade)
  [ OK ] 03 vRhyme        [vRhyme] OK (exit 0) -- rodou ate o fim. Bins formados: 0
  [ OK ] 03 VAMB          [VAMB] OK (exit 0) -- rodou ate o fim. Bins: 5
  [ OK ] 03 derep         5 vOTUs (representantes dereplicados)
  [ OK ] 04 PhaBox2       5 predicoes de lifestyle
  [AVISO] 05 VITAP         5 vOTUs mas 0 linhagens. Pode ser legitimo ...
  >>> VEREDITO: INTEGRO COM AVISOS — etapas completaram; revise os [AVISO].
```

---

## 12. CheckV no bin inteiro — N-linking (metodologia vRhyme)

**Problema:** um bin (vMAG) do vRhyme/VAMB pode ter vários scaffolds que são
fragmentos do **mesmo genoma**. Rodar o CheckV em cada scaffold isolado e pegar a
melhor completude **subestima** a qualidade do bin.

**Solução (recomendação do artigo do vRhyme, Kieft et al. 2022):** unir os
scaffolds do bin numa única sequência separada por um espaçador de **1000 N's**,
e rodar o CheckV nessa sequência — assim ele estima a completude do **bin inteiro**.
O run de N's introduz stop codons em todos os frames, impedindo que o prodigal-gv
(chamado internamente pelo CheckV) preveja genes atravessando as junções.

**Implementação:**

| Componente | Papel |
|---|---|
| [`scripts/link_bins.py`](workflow/scripts/link_bins.py) | une scaffolds de cada bin com `N`×1000 → 1 seq por bin (`linked_bins.fna`) |
| regra `link_bins` (`binning.smk`) | gera `03_bins_vrhyme_vamb/{sample}/linked_bins.fna` |
| regra `checkv_bins` (`checkv.smk`) | roda `checkv end_to_end` nos bins linkados → `02_checkv/{sample}_bins/quality_summary.tsv` |
| `filter_units.py` (modo metagenome) | usa a completude/qualidade **do bin** (chave `bin_id`), sem agregação por contig |

**Importante — só afeta bins multi-scaffold.** Bins de 1 scaffold dão resultado
idêntico (nada a unir). O espaçador é configurável em `config.yaml`:
`binning.link_spacer_n` (default 1000).

**Caminho sem reads (genoma de isolado):** continua usando o CheckV por contig
(não há binning), inalterado.

**Fluxo atualizado do bloco 03/02 (com reads):**
```
vRhyme + VAMB → collect_bins → link_bins (une com N's) → checkv_bins (CheckV no bin)
                                                              ↓
                                            filter_units (completude do bin inteiro)
                                                              ↓
                                                  derep (scaffolds reais) → vOTUs
```
> As sequências que seguem para derep / PhaBox2 / VITAP / abundância são os
> **scaffolds reais** (sem os N's). Os N's existem só para a estimativa do CheckV.

---

## 13. Anotação funcional: resistoma + mobiloma + enriquecimento

Bloco opcional (06_annotation / 07_annotation) que caracteriza **resistoma** (ARGs)
e **mobiloma** (elementos móveis) dos vOTUs e testa, por **permutação**, se esses
elementos estão **enriquecidos** na fração viral vs. o genoma completo.

### 13.1 Resistoma — ABRicate (TODOS os bancos + dedup com checkpoint)

[`abricate_resistome`](workflow/rules/annotation.smk) roda o ABRicate sobre os
vOTUs. Por padrão `abricate_dbs: "all"` → **TODOS os bancos embutidos** (resolvidos
em runtime via `abricate --list`; hoje 12: `ncbi`, `card`, `resfinder`, `argannot`,
`megares`, `vfdb`, **`victors`**, `ecoli_vf`, `upec_expec_vf`, `plasmidfinder`,
`ecoh`, **`bacmet2`** — biocidas/metais). Cada banco roda com `--threads` e tem a
saída crua guardada em `tools/abricate/<db>.tab` (proveniência/comparação).

> Os bancos **vêm embutidos e indexados** no pacote do ABRicate — não há download.
> O bloco 0 ([`setup_abricate_db`](workflow/rules/setup.smk)) valida (`abricate
> --check`) e grava um **manifesto** (`abricate_db_manifest.tsv`: nome, nº de
> sequências, tipo, data de cada banco) para reprodutibilidade. Para atualizar um
> banco à versão mais recente da fonte: `abricate-get_db --db <nome> --force`.

**Dedup por locus** ([`dedup_abricate.py`](workflow/scripts/dedup_abricate.py)):
hits da mesma sequência são agrupados por **sobreposição recíproca + union-find**
(cada intervalo cobre o outro em ≥ `dedup_overlap`). Isso funde "mesmo gene,
bancos diferentes" mas **não encadeia** genes adjacentes distintos nem funde um
gene aninhado num hit maior (corrige o encadeamento transitivo do clustering
guloso antigo). Mantém 1 linha por locus: o melhor hit (maior %id × %cov) +
todos os bancos/genes que concordaram (concordância = confiança).

Saídas: `resistome_abricate.tsv` (loci), `resistome_by_db.tsv` (contribuição de
cada banco: nº de hits brutos e em quantos loci participa) e
`resistome_dedup_checkpoint.json` — **checkpoint auto-verificável** (nº de bancos
rodados, bruto → loci, e `verification.status = PASS` só se não sobrou
sobreposição residual). O healthcheck (§11) lê esse checkpoint (`06 dedup-ARG`).

### 13.2 Mobiloma — múltiplas ferramentas consolidadas

Cada ferramenta emite **element calls** com coordenada no frame do vOTU
(`seq, start, end, element_type, tool, score`):

| Ferramenta | Sinal | Instalação |
|---|---|---|
| **geNomad** | provírus, AMR, conjugação (CONJScan), integrase | já no pipeline |
| **IntegronFinder** | integrons (complete / In0 / CALIN) | bioconda |
| **MobileElementFinder** | IS / transposons (MGEdb) | PyPI (env pip) |
| **CONJScan / MacSyFinder** | sistemas conjugativos (T4SS/relaxase) — **alternativa instalável ao ICEfinder** (módulo que o ICEfinder 2.0 adotou) | bioconda + `macsydata install --user CONJScan` (subconjuntos `CONJScan/Chromosome` + `CONJScan/Plasmids`) |
| **tRNAscan-SE** | tRNAs (marcador de sítio de integração) | bioconda |
| **ICEfinder** (opcional) | ICEs | manual (registro SJTU; não conda) |

> **ICEfinder vs CONJScan:** o ICEfinder é a referência de ICE, mas é web/registro,
> sem conda/pip — inviável num pipeline automático. O **CONJScan** é o motor de
> conjugação que o próprio ICEfinder 2.0 passou a usar e é instalável; por isso é
> a alternativa padrão (host-agnóstica). O ICEfinder fica como gancho manual opcional.

### 13.3 Dedup de elementos móveis — 2 CAMADAS + checkpoint

[`mge_dedup.py`](workflow/scripts/mge_dedup.py) — porque várias ferramentas
reportam o mesmo elemento físico no mesmo locus:

- **Camada 1 (intra-ferramenta):** funde calls redundantes dentro de cada tool
  (ex.: ABRicate multi-banco, CONJScan multi-modelo).
- **Camada 2 (entre-ferramentas):** funde calls de todas as tools no mesmo locus
  → 1 região não-redundante, anotada com **quais tools concordam** (concordância =
  confiança).

Critério **inteligente = sobreposição recíproca**: funde "mesmo elemento, tools
diferentes" (spans parecidos) mas **NÃO funde elementos aninhados** (um IS dentro
de um ICE permanece distinto). [`consolidate_mobilome.py`](workflow/scripts/consolidate_mobilome.py)
gera `mobile_element_regions.tsv` + `mobilome_dedup_checkpoint.json`.

**Checkpoint auto-verificável:** `verification.status = PASS` só se a camada 2
não deixou nenhuma sobreposição residual (`residual_overlaps = 0`). O healthcheck
(§11) lê esse checkpoint e acusa `FALHA` se a dedup não foi verificada.

### 13.4 Enriquecimento — permutação POR CATEGORIA

**Pergunta:** ARG/MGE estão sobre-representados na fração viral (prófagos/vOTUs)
vs. o genoma completo? Só faz sentido quando o genoma de entrada **tem fração
não-viral** (genoma bacteriano/MAG com prófagos); para fago puro o motor reporta
`NO_BACKGROUND`.

- [`build_gene_table.py`](workflow/scripts/build_gene_table.py): rotula TODOS os
  genes do geNomad (universo) como `is_viral` + **categorias**: `is_arg`,
  `is_conjugation`, `is_integrase`, `is_transposase`, `is_recombinase`,
  `is_integron`, `is_conjugative_system`, `is_IS`, `is_mge`. Para um teste de
  enriquecimento **válido** cada categoria precisa ser rotulada de forma idêntica
  em TODOS os genes do genoma (viral + não-viral). As fontes, **todas no genoma
  inteiro**, são:
  - **geNomad** (anota o genoma todo): `is_arg`, `is_conjugation`, `is_integrase`,
    `is_transposase`, `is_recombinase`;
  - **IntegronFinder** no assembly ([`integron_genome`](workflow/rules/enrichment.smk)) → `is_integron`;
  - **MobileElementFinder** no assembly ([`mefinder_genome`](workflow/rules/enrichment.smk)) → `is_IS`;
  - **CONJScan** no assembly ([`conjscan_genome`](workflow/rules/enrichment.smk), `--db-type unordered`) → `is_conjugative_system`;
  - `is_mge` = união de todas as anteriores.

  > As ferramentas mefinder/CONJScan rodam **2x**: nos vOTUs (mobiloma por fago,
  > §13.2) **e** no genoma inteiro (background do enriquecimento). Sem o segundo
  > uso, `is_IS`/`is_conjugative_system` ficariam sempre 0 (sem fração não-viral
  > para comparar) → permanentemente `UNDERPOWERED`. Por isso são duas execuções.

  **Qual fração viral entra (`is_viral`).** ARG/MGE são esperados em fagos
  **TEMPERADOS** (integram e trocam DNA com o hospedeiro), não em líticos nem em
  vírus eucarióticos. Por isso o default `enrichment_viral_set: temperate_phage`
  restringe `is_viral` a **fagos temperados**, com dois portões:
  - **Portão de fago (taxonomia geNomad):** mantém só vírus procarióticos
    (`phage_taxa_classes`, p.ex. *Caudoviricetes*) — **exclui vírus eucarióticos/RNA**
    (Orthornavirae, *Pararnavirae* = retrovírus tipo HIV, Megaviricetes…). geNomad
    detecta esses vírus eucarióticos, e o **PhaTYP assume que a entrada é fago** —
    então sem este portão ele rotularia um HIV como "temperate" sem sentido.
  - **Portão temperado:** provírus (integrados ⇒ temperados) **OU** contig de fago
    com **PhaTYP = temperate** (líticos são descartados).

  Modos: `temperate_phage` (default) · `provirus` (só integrados, mais estrito) ·
  `all_viral` (qualquer região viral; menos específico). O log do `build_gene_table`
  registra quantas regiões entraram e quantas foram filtradas (não-fago / líticas).
  **Trade-off:** restringir reduz o nº de genes virais → categorias raras podem
  virar `UNDERPOWERED`; em troca, o teste fica biologicamente específico e defensável.
- [`permutation_enrichment.py`](workflow/scripts/permutation_enrichment.py):
  **permutação estratificada por amostra** (preserva nº de genes virais → controla
  o confundidor de tamanho) + **hipergeométrica** como cruzamento. Auto-detecta
  TODAS as colunas `is_*` → 1 teste por categoria.
- **Guardrails (anti-conclusão-silenciosa):** `NO_BACKGROUND` (sem fração não-viral),
  `UNDERPOWERED` (poucos eventos — comum com ARG/integron raros), `OK`.

Saída: `enrichment_results.tsv` + `enrichment_checkpoint.json`.

### 13.5 Relatório consolidado

[`annotation_report.py`](workflow/scripts/annotation_report.py) gera:
- **`annotation_master.tsv`** — TSV geral por vOTU (mobiloma + nº loci ARG + tools
  concordantes por região).
- **`annotation_stats.log`** — log estatístico legível: enriquecimento por categoria
  (ratio, p permutação + hipergeométrica, status, **interpretação automática**
  ENRIQUECIDO/DEPLETADO/sem/UNDERPOWERED) + checkpoint da dedup + totais.

### 13.6 Configuração ([config.yaml](config/config.yaml) → `annotation`)

```yaml
annotation:
  abricate_dbs: "all"           # TODOS os bancos embutidos (ou lista p/ subconjunto)
  dedup_overlap: 0.5            # sobreposição recíproca p/ fundir loci/elementos
  run_integron_finder: true     # IntegronFinder nos vOTUs
  run_mobileelementfinder: true # MobileElementFinder nos vOTUs
  run_conjscan: true            # CONJScan/MacSyFinder nos vOTUs (alt. ao ICEfinder)
  run_trnascan: true
  run_icefinder: false          # gancho manual (icefinder_path)
  conjscan_models: "CONJScan/Chromosome CONJScan/Plasmids"  # subconjuntos do CONJScan 2.x
  run_integron_genome: true     # IntegronFinder no genoma inteiro -> is_integron
  run_is_genome: true           # MobileElementFinder no genoma inteiro -> is_IS
  run_conjscan_genome: true     # CONJScan no genoma inteiro -> is_conjugative_system
  enrichment_permutations: 10000
  enrichment_min_events: 5      # < isso -> UNDERPOWERED
  enrichment_viral_set: "temperate_phage"   # temperate_phage | provirus | all_viral
  phage_taxa_classes: "Caudoviricetes|Malgrandaviricetes|Faserviricetes|Leviviricetes"
```

---

## 14. Submissão ao ENA: pacote MIUViG dos uViGs

Bloco **modular e mode-aware** ([`ena.smk`](workflow/rules/ena.smk) →
[`build_ena_tables.py`](workflow/scripts/build_ena_tables.py)) que prepara,
**automaticamente**, o pacote de submissão dos uViGs (Uncultivated Virus Genomes)
ao **ENA** seguindo o padrão **MIUViG** (Roux et al. 2019, GSC) e a rota de
**assembly MAG** do Webin-CLI. Roda no fim (`stage_ena`; entra no `stage_summary`
quando `ena.enable: true`). Saídas em `0X_summary_final/ena_submission/`.

### 14.1 O que é auto-preenchido (do que o pipeline já produz)

| Campo MIUViG / Webin | Fonte no pipeline |
|---|---|
| `assembly_qual` | CheckV `miuvig_quality` (Finished / High-quality draft / Genome fragment(s)) |
| `completeness`, `compl_score`, `compl_appr` | CheckV (completeness + completeness_method) |
| `contamination` | CheckV |
| `source_uvig`, `detec_type` | `is_provirus` (mobiloma / CheckV / nome) → "provirus (UViG)" vs "independent sequence (UViG)" |
| `pred_genome_struc`, `topology` | geNomad/CheckV (DTR/ITR/Provirus) |
| `vir_ident_software`, `feat_pred`, `ref_db`, `sim_search_meth` | geNomad + versão (config) |
| `tax_class`, `vitap_lineage` | **VITAP** (taxonomia de saída) |
| `otu_class_appr` | derep (ANI/AF do config) |
| `number_contig`, `sequence_length` | vOTU FASTA |
| `COVERAGE` | abundância CoverM (metagenoma) **ou** `_cov_` no nome do contig (genoma) |
| `PLATFORM` | `read_type` das amostras (ILLUMINA / OXFORD_NANOPORE) |

### 14.2 O que VOCÊ preenche (só o submissor sabe) → `TODO_PREENCHER`

`collection_date`, `geographic_location`, `virus_enrich_appr`, `nucl_acid_ext`,
`assembly_software` (assembler a montante), `study_accession`, `sample_accession`,
`tax_id`. Preenchíveis via `config.yaml → ena:` (propagam pro pacote inteiro).

### 14.3 Arquivos gerados

- `ena_uvig_metadata.tsv` — tabela-mestre legível (todos os campos MIUViG por uViG).
- `ena_sample_registration.tsv` + `ena_samples.xml` — registro de 1 **SAMPLE por uViG**
  (cada MAG exige sample próprio), com checklist `ENA-CHECKLIST` e atributos MIUViG.
- `manifests/<alias>.manifest` — manifesto **Webin-CLI** (`ASSEMBLY_TYPE = Metagenome-Assembled Genome (MAG)`) por uViG.
- `fasta/<alias>.fasta.gz` — sequência do uViG com header limpo.
- `submit_all.sh` — loop que chama o `webin-cli` (validate → submit) em cada manifesto.
- `ENA_SUBMISSION_GUIDE.md` — passo-a-passo + softwares/versões citados.

O healthcheck (§11) confere consistência (uViGs = nº de manifestos) e conta os
campos `TODO` (item `09 ENA`). O exato **checklist** e **tax_id** são confirmados
no Webin (variam); o default é `ena.sample_checklist` (ver navegador de checklists
do ENA) — todo o resto é derivado automaticamente.

---

## 15. Resultado da validação — amostra MS13934

**Configuração do teste:**
- Amostra: `MS13934` (genoma completo de fago)
- Modo: `genome` (com reads paired-end Illumina)
- Plataforma: WSL2 Ubuntu 22.04, 8 GB RAM, 8 cores
- Data: 2026-06-24 (re-validação com vRhyme corrigido + N-linking + healthcheck)

**Resumo dos resultados (`per_sample_summary.tsv`):**

| Campo | Valor |
|-------|-------|
| vOTUs detectados | 5 |
| Lisogênico | 1 |
| Lítico | 1 |
| Lifestyle indefinido | 3 |
| Medium-quality (≥ 50% completude) | 1 |
| Low-quality | 4 |

**Tabela de vOTUs (`viral_bins_table.tsv`):**

| vOTU (contig) | Tamanho | CheckV | Completude | Lifestyle | Host (CHERRY) | Taxonomia PhaBox2 |
|---------------|---------|--------|-----------|-----------|--------------|-------------------|
| NODE_35 | 22.507 bp | Medium-quality | 57.3% | Lisogênico (1.0) | *Enterococcus termitis* (CRISPR) | Caudoviricetes |
| NODE_39 | 19.218 bp | Low-quality | 45.1% | Lítico (0.93) | — | Caudoviricetes |
| NODE_33 | 14.787 bp | Low-quality | 35.4% | — | — | — |
| NODE_48 | 13.114 bp | Low-quality | 33.5% | — | — | — |
| NODE_59 | 4.581 bp | Low-quality | 11.4% | — | — | — |

**VITAP:** sem atribuição taxonômica — fragmentos com < 60% de completude não
acumulam proteínas suficientes para superar os thresholds do grafo bipartido ICTV.
Resultado biologicamente esperado.

**Etapas concluídas com sucesso:**

```
✔ 01_genomad        identificação de contigs virais (5)
✔ 02_checkv_bins    CheckV no BIN inteiro (N-linking) — 5 vMAGs avaliados
✔ 03_binning        vRhyme (rodou até o fim, 0 bins legítimo) + VAMB (5 bins) + derep (5 vOTUs)
✔ 04_phabox2        lifestyle + host + taxonomia de alto nível
✔ 05_vitap          rodou até o fim (sem linhagem por fragmentação)
✔ 06_summary_final  tabelas finais integradas
✔ healthcheck       RUN_HEALTHCHECK.txt → "INTEGRO COM AVISOS" (0 FALHAS)
```

> **Nota sobre o vRhyme:** na primeira execução (2026-06-23) o vRhyme aparentava
> "0 bins"; a auditoria de 2026-06-24 revelou que ele **falhava em silêncio** (§10.8)
> e os 5 bins vinham só do VAMB. Após a correção, o vRhyme **roda até o fim** e
> confirma 0 bins de forma legítima (5 contigs / 1 amostra não têm sinal de
> co-cobertura; o próprio vRhyme avisa que opera melhor com 3+ amostras). Como o
> resultado de 0 bins se manteve, as tabelas finais permanecem idênticas e válidas.
> O `RUN_HEALTHCHECK.txt` agora atesta a integridade automaticamente.

---

## 16. Referências bibliográficas (ferramentas a citar)

Ao publicar resultados gerados com este pipeline, citar:

1. **geNomad:** Camargo AP, Roux S, Schulz F, et al. *Identification of mobile genetic elements with geNomad.* Nature Biotechnology, 2023. https://doi.org/10.1038/s41587-023-01953-y

2. **CheckV:** Nayfach S, Camargo AP, Schulz F, et al. *CheckV assesses the quality and completeness of metagenome-assembled viral genomes.* Nature Biotechnology, 2021. https://doi.org/10.1038/s41587-020-00774-7

3. **vRhyme:** Kieft K, Zhou Z, Anantharaman K. *vRhyme enables binning of viral genomes from metagenomes.* Nucleic Acids Research, 2022. https://doi.org/10.1093/nar/gkac341

4. **VAMB:** Nissen JN, Johansen J, Allesøe RL, et al. *Improved metagenome binning and assembly using deep variational autoencoders.* Nature Biotechnology, 2021. https://doi.org/10.1038/s41587-020-00777-4

5. **PhaBox2:** Shang J, Peng C, Guan J, Cai D, Wang D, Sun Y. *PhaBOX2: an enhanced web server for discovering and analyzing viral contigs in metagenomic data.* Nucleic Acids Research, 2026 (advance access). https://doi.org/10.1093/nar/gkag382

6. **VITAP:** Zheng K, et al. *VITAP: a high precision tool for DNA and RNA viral classification based on meta-omic data.* Nature Communications 16, 2226 (2025). https://doi.org/10.1038/s41467-025-57500-7

7. **Diamond:** Buchfink B, Reuter K, Drost HG. *Sensitive protein alignments at tree-of-life scale using DIAMOND.* Nature Methods, 2021. https://doi.org/10.1038/s41592-021-01101-x

8. **BLAST+:** Camacho C, Coulouris G, Avagyan V, et al. *BLAST+: architecture and applications.* BMC Bioinformatics, 2009. https://doi.org/10.1186/1471-2105-10-421

9. **Snakemake:** Mölder F, Jablonski KP, Letcher B, et al. *Sustainable data analysis with Snakemake.* F1000Research, 2021. https://doi.org/10.12688/f1000research.29032.2

10. **bwa-mem2:** Vasimuddin M, Misra S, Li H, Aluru S. *Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems.* IEEE IPDPS, 2019.

11. **samtools:** Li H, Handsaker B, Wysoker A, et al. *The Sequence Alignment/Map format and SAMtools.* Bioinformatics, 2009. https://doi.org/10.1093/bioinformatics/btp352
