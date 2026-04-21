## This script is the almost the same implementation as our notebook, just with class inheritance from the HMM.py class and a different example - the phmm_train_motif1.fasta file.
## It hasn't been altered so there's full integration with the HMM.py module - so complete 1 to 1 reassignment of the models we build with the attributes of the base class
## has not been performed. Additionally, no advanced analysis with training or classification has been done. This is just a simple example of how you can go about creating the pHMM topology
## and start to play around with integration of the original code.

from collections import defaultdict
from pprint import pprint

from HMM import HMM

class PHMM(HMM):
    """
    Class that contains all the methods for constructing a pHMM
    """

    def read_fasta(self, filename):

        sequences = []
        current_seq = ''

        with open(filename, mode='r', encoding='utf-8') as fh_in:
            for line in fh_in:
                line = line.rstrip("\n").rstrip('\r')
                # get the headers, reset current seq if one is found
                if line.startswith(">"):
                    if current_seq != '':
                        sequences.append(current_seq)
                    current_seq = ''
                else:
                    current_seq += line
        # get the sequence
        if current_seq:
            sequences.append(current_seq)

        return sequences

    def column_assignments(self, sequences):

        col_types = []
        num_seqs = len(sequences)
        seq_len = len(sequences[0])

        for i in range(seq_len):
            gaps = 0
            for j in range(num_seqs):
                if sequences[j][i] == '-':
                    gaps += 1

            gap_fraction = gaps / num_seqs
            if gap_fraction > 0.5:
                col_types.append("insertion")
            else:
                col_types.append("match")

        L = col_types.count("match")
        return col_types, L


    def sequence_states(self, sequences, col_types, L):

        seqs_states = []
        hidden_states = []

        for sequence in sequences:
            sequence_states = ["M_0"]
            counter = 1

            for j, char in enumerate(sequence):
                if char != '-' and col_types[j] == "match":
                    sequence_states.append(f"M_{counter}")
                    counter += 1
                elif char == '-' and col_types[j] == "match":
                    sequence_states.append(f"D_{counter}")
                    counter += 1
                elif char != '-' and col_types[j] == "insertion":
                    sequence_states.append(f"I_{counter-1}")
                elif char == '-' and col_types[j] == "insertion":
                    continue
                elif j == len(sequence) - 1:
                    sequence_states.append(f"M_{L + 1}")

            seqs_states.append(sequence_states)

        for seq in seqs_states:
            for state in seq:
                if state not in hidden_states:
                    hidden_states.append(state)

        self.hidden_states = hidden_states
        return seqs_states, hidden_states

    def initial_probabilities(self):
        return {"M_0": 1.0}

    def valid_transitions(self, L):
        transitions = defaultdict(list)

        for i in range(L + 1):
            if i < L:
                transitions[f"M_{i}"] = [f"M_{i + 1}", f"I_{i}", f"D_{i + 1}"]
                transitions[f"I_{i}"] = [f"M_{i + 1}", f"I_{i}", f"D_{i + 1}"]
                transitions[f"D_{i}"] = [f"M_{i + 1}", f"I_{i}", f"D_{i + 1}"]
            else:
                transitions[f"M_{i}"] = [f"M_{L + 1}"]
                transitions[f"I_{i}"] = [f"M_{L + 1}"]
                transitions[f"D_{i}"] = [f"M_{L + 1}"]

        return transitions

    def transition_probabilities(self, seqs_states, transitions, alpha):
        transition_probs = defaultdict(lambda: defaultdict(float))

        for state_path in seqs_states:
            for i in range(len(state_path) - 1):
                current_state = state_path[i]
                next_state = state_path[i + 1]
                transition_probs[current_state][next_state] += 1

        for state in transitions:
            for next_state in transitions[state]:
                transition_probs[state][next_state] += alpha

        for state in transitions:
            s = sum(transition_probs[state].values())
            for next_state in transitions[state]:
                transition_probs[state][next_state] = transition_probs[state][next_state] / s

        return transition_probs

    def emission_probabilities(self, col_types, seqs, L, b):

        emission_probs = defaultdict(lambda: defaultdict(int))

        # Match emissions
        idx = 1
        for i, state in enumerate(col_types):
            for seq in seqs:
                if state == "match":
                    if seq[i] != '-':
                        emission_probs[f"M_{idx}"][seq[i]] += 1
            if i != len(col_types) - 1:
                if col_types[i + 1] == "match":
                    idx += 1

        # Convert to probs
        for state in emission_probs.keys():
            s = sum(emission_probs[state].values())
            for emission in emission_probs[state].keys():
                emission_probs[state][emission] = (emission_probs[state][emission] + (1 * b)) / (s + (20 * b))


        # Insertion emissions
        background_probs = defaultdict(int)
        for seq in seqs:
            for emission in seq:
                if emission != '-':
                    background_probs[emission] += 1

        # Convert to probs
        total = sum(background_probs.values())
        for emission in background_probs.keys():
            background_probs[emission] = (background_probs[emission] + (1 * b)) / (total + (20 * b))

        # Assign to emission_probs
        for i in range(L + 1):
            emission_probs[f"I_{i}"] = background_probs.copy()

        return emission_probs


phmm = PHMM(hidden_states='MID', alphabet="ACDEFGHIKLMNPQRSTVWY")

seqs = phmm.read_fasta("phmm_train_motif1.fasta")
print("MSA Sequences:")
print(seqs)
states, l = phmm.column_assignments(seqs)
print("State Assignments by Position in MSA:")
print(states)
print("Total Number of Match States:")
print(l)
seqstates, hiddenstates = phmm.sequence_states(seqs, states, l)
print("State Assignments for Every Position of Every Sequence:")
print(seqstates)
print("Unique Hidden States:")
print(hiddenstates)
initial_probabilities = phmm.initial_probabilities()
print("Initial Probabilities:")
pprint(initial_probabilities)
transitions = phmm.valid_transitions(l)
transition_probs = phmm.transition_probabilities(seqstates, transitions, alpha=1e-3)
print("Valid Transitions:")
pprint(transitions)
print("Transition Probabilities:")
pprint(transition_probs)
emit_probs = phmm.emission_probabilities(states, seqs, l, b = 1)
print("Emission Probabilities:")
pprint(emit_probs)

