import sys
import numpy as np

def readFasta(path):
    # extract sequences from input file
    sequences = []
    with open(path, 'r') as file:
        sequence = ''
        for i in file:
            if i.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += i.strip()
        if sequence:
            sequences.append(sequence)
    return sequences[0]

def nucFreq(sequence):
    # calculate frequencies of each nucleotide base
    nuc_freq = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for s in sequence:
        for base in s:
            nuc_freq[base] += 1
    for base in nuc_freq:
        nuc_freq[base] /= len(sequence)
    return nuc_freq

def dinucFreq(sequence):
    # calculate frequencies of each dinucleotide
    totalLength = len(sequence) - 1
    dinuc_freq = {n1 + n2: 0 for n1 in 'ACGT' for n2 in 'ACGT'}
    for i in range(len(sequence) - 1):
        dinucleotide = sequence[i:i+2]
        dinuc_freq[dinucleotide] += 1
    for dinuc in dinuc_freq:
        if totalLength > 0:
            dinuc_freq[dinuc] /= totalLength
    return dinuc_freq

def expectedDinucFreq(nuc_freq):
    # calculate expected frequencies of each dinucleotide
    expected_freq = {}
    for n1 in 'ACGT':
        for n2 in 'ACGT':
            expected_freq[n1+n2] = nuc_freq[n1] * nuc_freq[n2]
    return expected_freq

def compareFreq(observed, expected):
    # compile observed and expected dinucleotide frequencies
    results = {}
    for d in observed:
        observed_freq = observed[d]
        expected_freq = expected[d]
        results[d] = {
            'Observed': observed_freq,
            'Expected': expected_freq
        }
    return results

def readIntervals(path):
    # stores all intervals from chrA.islands into an array
    intervals = []
    with open(path, 'r') as file:
        for i in file:
            start, end = map(int, i.strip().split())
            intervals.append((start, end))
    return intervals

def extractCpG(sequence, intervals):
    # extracts cpg regions and joins into 1 string
    joined = ''
    for start, end in intervals:
        joined += sequence[start:end-1]
    return joined

def extractNonCpG(sequence, intervals):
    # extracts non-cpg regions and joins into 1 string
    joined = ''
    current = 0
    for start, end in intervals:
        if current < start:
            joined += sequence[current:start-1]
        current = end
    if current < len(sequence):
        joined += sequence[current:]
    return joined

def markovModel(dinuc_freq):
    # builds 4 stage model based on frequencies
    transition = np.zeros((4,4))
    nucleotides = 'ACGT'
    index = {nuc: i for i, nuc in enumerate(nucleotides)}
    for d, freq in dinuc_freq.items():
        n1, n2 = d[0], d[1]
        i, j = index[n1], index[n2]
        transition[i,j] = freq
    for i in range(len(nucleotides)):
        total = sum(transition[i])
        if total > 0:
            transition[i] /= total
    return transition

def calcPotential(S, b, e, CpG, nonCpG):
    # calculate the cpg potential score over a window size
    nucleotides = 'ACGT'
    index = {nuc: i for i, nuc in enumerate(nucleotides)}
    prob_CpG = 1.0
    prob_nonCpG = 1.0
    for i in range(b, e - 1):
        n1, n2 = S[i], S[i + 1]
        j, k = index[n1], index[n2]
        prob_CpG *= CpG[j, k]
        prob_nonCpG *= nonCpG[j, k]
    if prob_nonCpG == 0:
        return float('inf')
    return np.log(prob_CpG / prob_nonCpG)

def predictIslands(S, CpG, nonCpG, windowSize=400, threshold=50):
    # predict islands based on cpg potential score
    predictions = []
    current = None
    for i in range(0, len(S) - windowSize + 1):
        b, e = i, i + windowSize
        score = calcPotential(S, b, e, CpG, nonCpG)
        if score > threshold:
            if current is None:
                current = b
        else:
            if current is not None:
                predictions.append((current, e - 1))
                current = None
    if current is not None:
        predictions.append((current, len(S) - 1))
    return predictions

def evalPredictions(predictions, actual):
    truePos = 0
    falsePos = 0
    falseNeg = 0
    evals = []
    # check for true positives (if (a, b) and (x, y) have significant overlap)
    for (x, y) in predictions:
        found = False
        for (a, b) in actual:
            # calculate overlap
            o = min(b, y) - max(a, x)
            if o > 0 and (4 * o) > (y - x):
                truePos += 1
                found = True
                evals.append((x, y, "TP"))
                break
        if not found:
            falsePos += 1
            evals.append((x, y, "FP"))
    # check for false negatives
    for (a, b) in actual:
        found = False
        for (x, y) in predictions:
            # calculate overlap
            o = min(b, y) - max(a, x)
            if o > 0 and (4 * o) > (y - x):
                found = True
                break
        if not found:
            falseNeg += 1
            evals.append((a, b, "FN"))
    return truePos, falsePos, falseNeg, evals

def main():
    # extract sequences (and intervals) from files
    fasta = sys.argv[1]
    flags = sys.argv[2:]
    intervalsFile = sys.argv[-1] if '-hmm' in flags else None
    sequence = readFasta(fasta)
    
    # specify options
    freq = '-f' in flags
    hmm = '-hmm' in flags
    
    if freq:
        # perform frequency calculations (observed and expected)
        nuc_freq = nucFreq(sequence)
        dinuc_freq = dinucFreq(sequence)
        expected = expectedDinucFreq(nuc_freq)
        # compare between observed and expected frequencies
        comparison = compareFreq(dinuc_freq, expected)
        # identify significantly different dinucleotides
        significant = {k: v for k, v in comparison.items() if 
                   abs(v['Expected'] - v['Observed']) > 0.02}
    
        # write calculations to stdout
        print("Frequency Comparison:")
        print(f"{'Dinucleotide':<12}{'Observed':<10}{'Expected':<10}")
        for dinuc, values in comparison.items():
            print(f"{dinuc:<12}{values['Observed']:<10.4f}{values['Expected']:<10.4f}")
        print("Significant Dinucleotides:")
        print(" ".join(significant.keys()))
    
    if hmm:
        # read intervals and extract CpG and non-CpG regions
        intervals = readIntervals(intervalsFile)
        CpG = extractCpG(sequence, intervals)
        nonCpG = extractNonCpG(sequence, intervals)
        # calculate frequences for CpG and non-CpG
        CpG_dinucFreq = dinucFreq(CpG)
        nonCpG_dinucFreq = dinucFreq(nonCpG)
        # make markov models and predictions for CpG and non-CpG
        CpG_model = markovModel(CpG_dinucFreq)
        nonCpG_model = markovModel(nonCpG_dinucFreq)
        predictions = predictIslands(sequence, CpG_model, nonCpG_model)
        
        # output tp/fp/fn summary and predictions
        truePos, falsePos, falseNeg, evals = evalPredictions(predictions, intervals)
        print(f"True Positives: {truePos}")
        print(f"False Positives: {falsePos}")
        print(f"False Negatives: {falseNeg}")
        with open('results.txt', 'w') as out:
            out.write("Start\tEnd\tClassification\n")
            for start, end, eval in evals:
                out.write(f"{start}\t{end}\t{eval}\n")

if __name__ == "__main__":
    main()