#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$(pwd)}"
LOGDIR="${LOGDIR:-$HOME/log}"
JOB_NAME="${JOB_NAME:-segtrace-angio-genus}"
THREADS="${THREADS:-128}"
MEM="${MEM:-480gb}"
WALLTIME="${WALLTIME:-96:00:00}"

mkdir -p "$LOGDIR"
cd "$WORKDIR"
mkdir -p selected results selected/clean_fasta

DATA_DIR="./eukaryotic_data/ncbi_dataset/data" \
SELECTED_DIR="./selected" \
python3 - <<'PY'
import base64
import gzip
import json
import os
import re
import shutil
import sys
import tempfile
import time
import urllib.request
import zipfile
import zlib
from pathlib import Path

# This is the 466-species list from AngioWGD's "Original data" table, embedded
# here so this submit script has no runtime dependency on a separate list file.
SPECIES_B85 = "c-n1R+jiWz5q#%Y@CWh{qmkBH$FdZU;^g?PNkSwf5a0sDtf;R~RfFV=eDaXBfZ?Ulm#XS!4r3FdaZc(v+I9$2=%o6-jME?ndxMFi#5`WM+TcQzS*p7n7jdS3F)_SXH&K%DR_fFa(KII7>OREI<gq?y*Gn|j$+Wz_HPPVT>N(okw)lO{v9}%GaWK7gVM_HZ8IQ5{a*dcY+cU0t9;^@30+VX%XIN)C%UjRJ3wBa=5^Y-7NpO?Uwsu%>GF({?OB3~NpLRQ&W-f!T4(;HD2gi1robjFUsk#X~oQ&9bx2<G{mwd_6Wo-IB!@8H?tW=L~iBdb+)<2}{aATjXcsvJpJ;voj>X={b6nea`tXEUZqIUuhMmugBQ}*plr++ghCOmKA&vf!H`ysUK^2Tjv5(l&Bj0c-q@$-;V?eYx&glH%HlT&q)wjCJYI=sBt)_M*0Da3xjUSl=dIQ%@|u6A%i>dv%wvh>L>q{Apns&0`ve1UQFC9LQ}(siotS_Z7#*~Xmy(Tj5yxj;LlK3fC>CT?pxi-UWsC*!?MX2#>-jEoeho^UFr0&eLAKfpC%GBe*e4`CFfTJ<X|^v<hH2p1l+u1A?Q$E&+2eGCW|a?HG3FTK?`!<451yQ7^=s)Y^|r<-t0+PzHA;gQpflL<`!oC80G#iNwcghHvCXwJC4MPh#bVakZ{Z+o-x8E-qHT;Aa>CgrS#mz~rG#!afP$Oo8XUS4v7E$h<y7i1^a2QOX{N_7nU`9C?nG-AE_+d0Ih`X>GzHe+W_65ax`1}hk}=K(78&LaY*Mkx^+@ThUy+G6sM6HX|uv2e4OvKIV;2h6nE$q2)eMcv@3tLtvXjKPP<Qhzo64P~o$9|K-e`-}_enL8&2kwNKrQ$llyfRMU{Pud|PCIN`^UU$;Yi!1jaPPhvQfy-<2S2m>rj?s+5Sk^L4sA$V@z`$Gl9Yez$^7sd5K1|35VbY|*_z(P#d311^*{zxOzqMeAJP9EYME~MBTDtGTU(y-|MQf0yAmfsd7ckq<i01;lxc0*=PP3dP6%cJ8WgWwueW%v$d!p7qg;~bXn>lzoeNW9_I0Br3w2v_7vdl;&V#aY9Ll;^+$K&6`2PNAnz;Rbk{7gqC`T>P7=gxZK!enFU49H>%z_|fPbG5hQF=1*1{2~hg9{b;bp6Wlr!Ar7N7Q_C8G^f9@G24LGRfjRe2`+&Rz_9u)_*ZGGU+|yl44FY{90k|T`hWwRlF1sx5Fh>qajaW1;=J%Xobm+Vr03Qh(e!Cd-Gm&yfm5va@aBlnz5%KS8^u(|!6H*v>`Zc%!`}d0I5yMAgTWjiQ4{dfKCus(YpMt7J)r@J8jyV!X9j`;bb5tb$uZO1^?J>8I0e|NXMM3FIExYgwA}_jhsZJ;NXe^_20?S1m7ji@aPFC#uyclU7dRKpXW$W*_AyKZ!ix<|0o%?|0k-@ZY{&5Yl~F3angE(~My|HP`;t&5PQ&8utcdLpP~s1h=x(xwv=sqUOD6dz(5EC|_V$G~SBE##SXCtos!Vs8<P5K+FnT&$Ibq!b5A>gQBBCpxK51dzLQA98&@sfNV~N&t7y)dGuFnRt1Idvla2sj?E(G=BUg>1CCA0h0kgM%i@zgwsCTf(@XaRi|7kxMIZ4)kF_H2;eEs+#?k=C?$j)>%*Ey98*x8LM<g-tERve1eNOGcKi@%au+)&f|M0>fWI96$de4dM7|M`n_Y`g2&4dICp7X<XJ*1mqGI{)BlE)QXegrOXilu&N&OK3}i*Oy76#-d1V)6ddCY$>Ke$*6urQ0iRLNv{#qSq($YLTihNCE}#gh@}3|FREE8S##HAW$*eOuVVpLsU$-vg{pSQCB-&ZWiU!NTG+<~Wta?Pi%T+EOf#bLs`loTlXWJHF!vK~L%`XCJ>oo-uQErVnjsVPU-&h&U5xn4W!V9B9CEhVw<#LqMR+Ts^$pK6m5Z;N=uswl}P%={&08`->6BlT(zz0Dyo^V2=E+G9{>gu9hpeO<MBe6#sXFc++oq4_<FQm}=5t7ovZ+vwUcP0k!diOEXj4SoNZ;bLvOd&-qBzIK42hq0VO~T7#Ye$VgaY%(55SI<cON7y@k0}CMiPL95hGGPnY)4*I9CKsl%WF0xz`TJ)v2Ky6XBCVTpOYY(2F*gh0LxIP)QC;l#B=zPZ>sZ<V#F$!p(EC6dV;18e>0{iBBU88;Q<FB|MYn-NBfFL(vPIGQOk{ha4)RiHJ(s9#3I&@3Nicw*CPo~R!Z&4+_0W<qiAH9hEt$ru`N~a&f6!GWYj*VQ!o6RDlC_vL#rXVz8LSUvRZ~J08n(DGp|f-FtK0FJNVf$0oT}FNmuxv8Hk3)R)0hk6|;h6nih91GT;&DS(}W^MDJIf78nMvhel8FRg+(}5K_IZ{)GBeWdN**s#%ffOpzH?2C5<>9b&7J-9d=q5+WvxC3cjsj>e71jg7JbLbnTUW;q3k5J{o<O4Z#PstCPYL5c0on>IkBZMDkiZo}+GUn)>COoO{fzR)u9(ruiJ%CY(eJczU!m(s>C3ZxhWYBJRqXGr7ZN@!TGg8UmS%I)69J>yjG1YT&Dw|^w*hvnmH6ItAcQXEHWzetZzmcq)b9nj9w;6U1Y>uF0aAlu=8fF~$(Eo|jXCb5_Its+;ymA$Q)wly*o+|{ou2(XDf?lsL1s>~p<w9R)Wvl~{_RDpqfYS--U4a<1tGh#0cXEj{w<B3Osm1n>sYyXz5Lt@MW8(LaLb-MR#-<OoG9$A*xymtiF&*(_m3tW1QTRGcKx%wPX*Ks9#ru5&|C5j8g{eJ1-uU;G1HWm?tI<h?P65^YZIm;fclF>l|L=ElJV6XLyQ~d<HS+6~TpiHAx^1IrF%mDwB#(sd|caHDbwHIT4{|TBS;q^yn0Q5Ln<Y0y>3e}@0GtriGzHD8@+K<4g2F&i!yYu<#C-N=I)tMqf_8kBsNQy~l$n$J0rUl4bgM|>-pLS;#wlQ6J2e`1!Eo6Y#>1u*2sPlkJ9E;RQ+j$Mqz#goOR+x=J`vfm-v4%>BxalAk9hJ{GZf>YfSxhG?hJa8iTG0>eeA&i~1Z47E_STWkZ+jE96`lzm`1}*ZDLLnQOGQVCEnBbcT@sg7C#UMmG9@f5o1ZHnhwAm_YWng-@^He`1lhwC_c$R(UU;6iaYQnU^ZG7*3g9!P1G?F(kPe)l%^WT05Ulz``gQbGS{<cE*ecZV3#bs<)G74s5l19}Q@Hk{6jxcFtuH_Slh?ydoZ@veIxBSJ!D?Pwh#7brTCp6|bd-Ir!8`9^DiTQ_1RgDw^(3ZrkiZgK-9dyE*4yj)n*edCVNH?L3S_4MQOLe;)niv{1fPP6+4c+u6`e+*1c6|g4^-sSVtfEnp)G<=Wi_=0fM5a$13G{NzKWmI^0D;HNujDOPUX#jxb`>J{DW9o;{cQGwiO9QS(A;}&Yzwh%o{j8Xy54pZH<QiCYXWQM)i4S(&wTnT%MHfj^G+xnl|QmLQbZs3>%?T6e~CfH8@=$%*>q3{IVZ*Vv4@z*9^Da+i7CEY|!_}xav{PPtbQNF<wol!1R47Z9Z&8z<u)!g~sFcgn70Lz%^AP5RRYmCUJ<#L0|PFH<XjWzE17q&$@+pDVJ>w9wK+A0C%Nf_Za7tKtFS4D`3}khn6N4Ww+3yjbW3$)>zqzh8YW?*0Gww=%ImzUaqL*NER@h!c6S$;Q<EoDdr+{*m>+raJ;E*aXa3P^chQ^egl>7k%_XW!jXfoEnW=Zx<)|+a5vgzomse*OZjtbEc=nzFRQqftL+iGYb20Fu}pgoB+o$4sBc#&uV0{3)7BcUdmYu)Svq9}(%;M}Z{0gMPRZ{8(NHutC0y5=zX44gs((S>g)l;0y~@6T862Lm<Y4`d^FM9x%XHIStw9Oek(X7<4GIS(#-0UThcxJPP#VNMcFW*-VBom-_S!eTawvlu$G@NkbVR`gCPyL?qi5l*@DV$+n3Z0uL&LC7oB^(4^1%>9`mhzxzD5*$Ti6NpZI#s&yoe<w9mgxX`CVb(p5Q%Y^KbYsw7KI@97orsiRTiApO##<6|GI1;!#CsA|V%oktF~4zO0>{cGIolSc}4kHeF=d0Vy}xX<oU<>(J_3p6qW0T=ZE>XeGz>)aBbom*ZAdw>WKSav(8cLHh$|XN%{sz|<%=Pc#tQLESY4fSKTjR%v=#T3_VG9cmhkJeBMy$3LM<;M{?t0?heoZgeAgD8pttp)Zq4xAfxH0e7ilv<IkFbcXDo*IU2l9)&|R=-D!|$5?kFhfDhktz4Ie9JR=omP&yS1G8HW(pG(}dueU8s5l<hPhBxLv5bP4?1MoXY*G+b6_>y0cpR1Wm=si^8Ma!xAr>iqQgPvyl(l0kV{!Ua2fu8h6v6(LDK6OSH#uHbYQ(nE`E8jy;dwKQ32g$LQ}%8Ewb@c9C-Yv)5#$YmP<_{gGkvvhwE72V9JH?ki&~~lvD_o@{3@aC`G&yW4q>hhzfusvtDhvow)Db;!LPy;9`NKGfrix&4h|IH$1#(V?=s)hln+tE61Aa*>DrfoxO&;!?0pxCf>iO$`D;`dDHa}wP(CF>0_TmjP6+M$|JmX97iRRDmrp@yb(8#d-zZ-<3wb{n1dQJ;qfS}a!J-47r39dW!GxSBS^<w*y>Qgh)w+FQXxhuiMl{$#`&e|s5G{DEZ*45^-R74{4rv4Yb#-Z-%DQ!CEYz#>?tuzf1<R;X=X$TnY-9*A&w(SJ{ZS80T~_?wF6!4!Wp||93Gaw&#WXgSCIgVQ;I^hTGP<6U$z1?3<3wF)SPq&ZpC;A+EnhXXMM_M80;AP;IIDaQ)K<SG$lls<Q%#Wl$O+E4@TjZDosg+0EAsIhx@!ca`jfy~THiisLGfcVxx1z=ysW(7$(QBcX#9Z6QM%bU`z~Jf>L>V`=&btlF)sWTqn`lDgH*x#_2$IS+WyTdexNBjvjU%SLWd~)zU$bxAZBNE=3Smuz9n%fLoU-350sMgdLc<@6QjLsl+AA{XIU(mqFJXt>=0%q&ux@O6=nI@6jecrV~eHFD8tEV=||E7G$@pP8Kv5GUi(0M`9RL0?FNiG5Z5Pm<trLLXwRRw^m9NNqbyo6Y$&ZZQ>6^Q8EB>JuVjE#&dk9(gN=z_f(~eNZSZuA@X<O9K|4B{&FeZwF`m<VP)~9z`W~51(dq{col@Ui3&HDQ1Y0DAWx^wf{bhEwwMM5|KZjA&qMtFG<eKswLw&tJ<}~OhH)ke}rJeEyRjtaGY5PAMAQeJWX61-rsBS?{q5iN9k4Y-OH;f(M=j=@4WaJ*#`8v8fc&DSZzX@X_z%^qz;6ni{IzTCIgO4s-4bV_$SLSLHDsq*s_ZCS}CN0`x(1-)VI|ew5ACFhG$M*jK^9|0P"
SPECIES = zlib.decompress(base64.b85decode(SPECIES_B85)).decode().splitlines()
ALIASES = {"Pisum sativm": "Pisum sativum"}
LEVEL_SCORE = {"complete genome": 3, "chromosome": 2, "scaffold": 1, "contig": 0}
ORGANELLE = re.compile(r"^(chr)?(m|mt|mit|mito|mitochondria|mitochondrion|pt|cp|pltd|plastid)$", re.I)
ORGANELLE_WORDS = ("chloroplast", "mitochondrion", "mitochondrial", "chloroplastic", "plastid", "plastome", "chondriome")

data_dir = Path(os.environ["DATA_DIR"])
selected_dir = Path(os.environ["SELECTED_DIR"])
report_path = data_dir / "assembly_data_report.jsonl"
clean_dir = selected_dir / "clean_fasta"
out_path = selected_dir / "angio_wgd_genomes.files"
summary_path = selected_dir / "angio_wgd_genomes.tsv"
clean_dir.mkdir(parents=True, exist_ok=True)

def key(record):
  return (LEVEL_SCORE.get(record["level"].lower(), -1), record["is_ref"], record["busco"], record["scaffold_n50"], record["contig_n50"], record["size"])

records_by_name = {}
with open(report_path) as report:
  for line in report:
    data = json.loads(line)
    stats = data.get("assemblyStats", {})
    accession = data.get("accession") or data.get("currentAccession")
    name = data.get("organism", {}).get("organismName", "")
    size = int(stats.get("totalSequenceLength") or stats.get("totalUngappedLength") or 0)
    if not accession or not name or not size:
      continue
    busco = stats.get("busco", {})
    record = {
      "accession": accession, "name": name, "size": size,
      "level": data.get("assemblyInfo", {}).get("assemblyLevel", ""),
      "is_ref": int(data.get("assemblyInfo", {}).get("assemblyCategory", "").lower() in ("reference genome", "representative genome")),
      "busco": float(busco.get("complete") or busco.get("buscoScore") or 0),
      "scaffold_n50": int(stats.get("scaffoldN50") or 0), "contig_n50": int(stats.get("contigN50") or 0),
    }
    records_by_name.setdefault(ALIASES.get(name, name).casefold(), []).append(record)

selected = []
missing_species = []
for species in SPECIES:
  candidates = records_by_name.get(ALIASES.get(species, species).casefold(), [])
  if not candidates:
    missing_species.append(species)
  else:
    selected.append(max(candidates, key=key))
if missing_species:
  raise SystemExit("species absent from NCBI report: " + ", ".join(missing_species[:20]))
if len({record["accession"] for record in selected}) != len(selected):
  raise SystemExit("multiple AngioWGD species resolved to the same NCBI accession")

def genome_paths():
  return {path.parent.name: path for path in data_dir.glob("*/*_genomic.fna*")}

def download(accession):
  url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA"
  print(f"[download] {accession}", file=sys.stderr)
  with tempfile.TemporaryDirectory() as temporary:
    archive = Path(temporary) / f"{accession}.zip"
    for attempt in range(3):
      try:
        urllib.request.urlretrieve(url, archive)
        break
      except Exception:
        if attempt == 2:
          raise
        time.sleep(5 * (attempt + 1))
    with zipfile.ZipFile(archive) as zipped:
      members = [member for member in zipped.namelist() if member.endswith("_genomic.fna") or member.endswith("_genomic.fna.gz")]
      if not members:
        raise RuntimeError(f"no genomic FASTA in NCBI download: {accession}")
      destination = data_dir / accession
      destination.mkdir(parents=True, exist_ok=True)
      for member in members:
        with zipped.open(member) as source, open(destination / Path(member).name, "wb") as target:
          shutil.copyfileobj(source, target)

paths = genome_paths()
for record in selected:
  if record["accession"] not in paths:
    download(record["accession"])
paths = genome_paths()
missing_fastas = [record["accession"] for record in selected if record["accession"] not in paths]
if missing_fastas:
  raise SystemExit("missing NCBI genomic FASTAs: " + ", ".join(missing_fastas))

def is_organelle(header):
  text = header[1:].strip()
  return ORGANELLE.fullmatch(text.split()[0]) is not None or any(word in text.lower() for word in ORGANELLE_WORDS)

def filter_organelle(source, destination):
  opener = gzip.open if source.suffix == ".gz" else open
  temporary = destination.with_suffix(destination.suffix + ".partial")
  kept = discarded = 0
  write = False
  with opener(source, "rt", encoding="utf-8", errors="replace") as input_handle, open(temporary, "w") as output_handle:
    for line in input_handle:
      if line.startswith(">"):
        write = not is_organelle(line)
        kept += write
        discarded += not write
      if write:
        output_handle.write(line)
  os.replace(temporary, destination)
  return kept, discarded

with open(out_path, "w") as output, open(summary_path, "w") as summary:
  summary.write("species\taccession\tsize_bp\tassembly_level\traw_fasta\tused_fasta\tdiscarded_organelle_count\n")
  for record in sorted(selected, key=lambda item: item["name"]):
    raw = paths[record["accession"]]
    clean = clean_dir / f"{record['accession']}_{raw.name.removesuffix('.gz')}.nuclear.fna"
    if clean.exists():
      discarded = 0
    else:
      _, discarded = filter_organelle(raw, clean)
    output.write(f"{clean}\n")
    summary.write(f"{record['name']}\t{record['accession']}\t{record['size']}\t{record['level']}\t{raw}\t{clean}\t{discarded}\n")
print(f"selected_species={len(selected)}")
print(f"written_fastas={len(selected)}")
PY

qsub -N "$JOB_NAME" \
  -l select=1:ncpus=${THREADS} \
  -l walltime=${WALLTIME} \
  -v WORKDIR="${WORKDIR}",THREADS="${THREADS}" \
  -j oe \
  -o "$LOGDIR/${JOB_NAME}.log" <<'PBS'
#!/usr/bin/env bash
set -euo pipefail

cd "$WORKDIR"
mkdir -p results

mapfile -t FASTAS < ./selected/angio_wgd_genomes.files
printf '[segtrace] input FASTA count: %d\n' "${#FASTAS[@]}"

./time -v ./segtrace \
  -p "$THREADS" \
  -c 1 \
  -o ./results/ANGIOSPERM_GENUS \
  "${FASTAS[@]}"
PBS
