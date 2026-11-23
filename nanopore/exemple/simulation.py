import gzip, random, re
random.seed(7)
# --- références : 400 pb de vraies séquences 16S (jeu d'exemple dada2)
seqs = []
with gzip.open("example_train_set.fa.gz  # fichier extdata du package dada2", "rt") as f:
    s = None
    for l in f:
        l = l.strip()
        if l.startswith(">"):
            if s: seqs.append(s)
            s = ""
        else: s += l.upper()
    seqs.append(s)
seqs = [s for s in seqs if len(s) > 1000 and set(s) <= set("ACGT")]
base = [s[500:900] for s in random.sample(seqs, 4)]
# E : variante proche de A (8 différences, ~2 %)
E = list(base[0])
for i in random.sample(range(20, 380), 8):
    E[i] = random.choice([b for b in "ACGT" if b != E[i]])
refs = {"spA": base[0], "spB": base[1], "spC": base[2], "spD": base[3], "spE_procheA": "".join(E)}
F, R = "AGGTGAAGTAAAAGGTTCATACTTAAA", "CCTTCTAATTTACCAACAACTG"
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
ADAPT = "TTTCTGTTGGTGCTGATATTGC"

def mutate(s):
    out, qual = [], []
    i = 0
    while i < len(s):
        b = s[i]
        # longueur de l'homopolymère qui commence ici
        j = i
        while j < len(s) and s[j] == b: j += 1
        hp = j - i
        for k in range(hp):
            r = random.random()
            pdel = 0.006 + (0.03 if hp >= 4 else 0)
            if r < pdel: continue                                   # délétion
            if r < pdel + 0.006:                                    # substitution
                out.append(random.choice([x for x in "ACGT" if x != b])); qual.append(random.gauss(9, 3)); continue
            out.append(b); qual.append(random.gauss(24, 6) if hp < 4 else random.gauss(16, 5))
            if random.random() < 0.005 + (0.01 if hp >= 4 else 0):  # insertion
                out.append(b if hp >= 4 else random.choice("ACGT")); qual.append(random.gauss(8, 3))
        i = j
    q = "".join(chr(int(max(3, min(40, x))) + 33) for x in qual)
    return "".join(out), q

def read_from(seq):
    amp = ADAPT[random.randint(0, 10):] + F + seq + rc(R) + rc(ADAPT)[:random.randint(5, 22)]
    if random.random() < 0.5: amp = rc(amp)
    return mutate(amp)

def sample(path, comp_, n):
    names = list(comp_.keys()); w = list(comp_.values())
    with gzip.open(path, "wt") as f:
        for k in range(n):
            r = random.random()
            if r < 0.03: s, q = mutate(rnd(random.randint(100, 900)))                       # bruit
            elif r < 0.05:                                                                    # chimère
                a, b = random.sample([refs[x] for x in names if x in refs], 2) if len(names) > 1 else (refs[names[0]], refs["spB"])
                s, q = read_from(a[:200] + b[200:])
            else: s, q = read_from(refs[random.choices(names, w)[0]])
            f.write(f"@read{k}\n{s}\n+\n{q}\n")

sample("reads/barcode01.fastq.gz", {"spA": 45, "spB": 25, "spC": 15, "spE_procheA": 12, "spD": 3}, 3000)
sample("reads/barcode02.fastq.gz", {"spB": 50, "spC": 30, "spD": 20}, 2000)
sample("reads/barcode03.fastq.gz", {"spA": 100}, 800)
with open("mock_attendu.fa", "w") as f:
    for k, v in refs.items(): f.write(f">{k}\n{v}\n")
print({k: len(v) for k, v in refs.items()})
