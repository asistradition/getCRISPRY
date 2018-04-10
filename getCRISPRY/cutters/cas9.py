from __future__ import print_function
import numbers

from . import Cutter
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein


class cas9(Cutter):
    # Mismatch, insertion & deletion data from Doench et al 2016

    name = "spCas9"
    seq = Seq(
        "MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGALLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHRLEESFLVEED"
        "KKHERHPIFGNIVDEVAYHEKYPTIYHLRKKLVDSTDKADLRLIYLALAHMIKFRGHFLIEGDLNPDNSDVDKLFIQLVQTYNQLFEENPINASGVDAKAILSARLSKSR"
        "RLENLIAQLPGEKKNGLFGNLIALSLGLTPNFKSNFDLAEDAKLQLSKDTYDDDLDNLLAQIGDQYADLFLAAKNLSDAILLSDILRVNTEITKAPLSASMIKRYDEHHQ"
        "DLTLLKALVRQQLPEKYKEIFFDQSKNGYAGYIDGGASQEEFYKFIKPILEKMDGTEELLVKLNREDLLRKQRTFDNGSIPHQIHLGELHAILRRQEDFYPFLKDNREKI"
        "EKILTFRIPYYVGPLARGNSRFAWMTRKSEETITPWNFEEVVDKGASAQSFIERMTNFDKNLPNEKVLPKHSLLYEYFTVYNELTKVKYVTEGMRKPAFLSGEQKKAIVD"
        "LLFKTNRKVTVKQLKEDYFKKIECFDSVEISGVEDRFNASLGTYHDLLKIIKDKDFLDNEENEDILEDIVLTLTLFEDREMIEERLKTYAHLFDDKVMKQLKRRRYTGWG"
        "RLSRKLINGIRDKQSGKTILDFLKSDGFANRNFMQLIHDDSLTFKEDIQKAQVSGQGDSLHEHIANLAGSPAIKKGILQTVKVVDELVKVMGRHKPENIVIEMARENQTT"
        "QKGQKNSRERMKRIEEGIKELGSQILKEHPVENTQLQNEKLYLYYLQNGRDMYVDQELDINRLSDYDVDHIVPQSFLKDDSIDNKVLTRSDKNRGKSDNVPSEEVVKKMK"
        "NYWRQLLNAKLITQRKFDNLTKAERGGLSELDKAGFIKRQLVETRQITKHVAQILDSRMNTKYDENDKLIREVKVITLKSKLVSDFRKDFQFYKVREINNYHHAHDAYLN"
        "AVVGTALIKKYPKLESEFVYGDYKVYDVRKMIAKSEQEIGKATAKYFFYSNIMNFFKTEITLANGEIRKRPLIETNGETGEIVWDKGRDFATVRKVLSMPQVNIVKKTEV"
        "QTGGFSKESILPKRNSDKLIARKKDWDPKKYGGFDSPTVAYSVLVVAKVEKGKSKKLKSVKELLGITIMERSSFEKNPIDFLEAKGYKEVKKDLIIKLPKYSLFELENGR"
        "KRMLASAGELQKGNELALPSKYVNFLYLASHYEKLKGSPEDNEQKQLFVEQHKHYLDEIIEQISEFSKRVILADANLDKVLSAYNKHRDKPIREQAENIIHLFTLTNLGA"
        "PAAFKYFDTTIDRKRYTSTKEVLDATLIHQSITGLYETRIDLSQLGGD",
        alphabet=generic_protein)

    PAM = "NGG"
    PAM_3prime = True
    PAM_grid = {_key: {
        "A": {"G": 0.259259259},
        "T": {"G": 0.038961039},
        "G": {"A": 0.069444444, "T": 0.016129032, "G": 1, "C": 0.022222222},
        "C": {"G": 0.107142857}
    }
        for _key in ["A", "T", "G", "C"]}

    spacer_length = 20
    spacer_grid = {
        1: {"A": {"T": 1, "G": 1, "C": 0.857142857},
            "C": {"T": 1, "G": 0.913043478, "A": 1},
            "G": {"T": 1, "C": 0.714285714, "A": 0.9},
            "U": {"G": 0.956521739, "C": 0.857142857, "A": 1}},
        2: {"-": {"T": 0.692307692, "G": 0.714285714, "C": 0.96, "A": 0.727272727},
            "A": {"-": 0.969230769, "T": 0.727272727, "G": 0.8, "C": 0.785714286},
            "C": {"-": 0.96875, "T": 0.909090909, "G": 0.695652174, "A": 0.727272727},
            "G": {"-": 0.921875, "T": 0.636363636, "C": 0.692307692, "A": 0.846153846},
            "U": {"-": 0.953125, "G": 0.84, "C": 0.857142857, "A": 0.846153846}},
        3: {"-": {"T": 0.4375, "G": 0.461538462, "C": 0.666666667, "A": 0.6},
            "A": {"-": 0.888888889, "T": 0.705882353, "G": 0.611111111, "C": 0.428571429},
            "C": {"-": 0.933333333, "T": 0.6875, "G": 0.5, "A": 0.866666667},
            "G": {"-": 0.822580645, "T": 0.5, "C": 0.384615385, "A": 0.75},
            "U": {"-": 0.921875, "G": 0.5, "C": 0.428571429, "A": 0.714285714}},
        4: {"-": {"T": 0.15, "G": 0.176470588, "C": 0.375, "A": 0.1},
            "A": {"-": 0.580645161, "T": 0.636363636, "G": 0.625, "C": 0.352941176},
            "C": {"-": 0.682539683, "T": 0.8, "G": 0.5, "A": 0.842105263},
            "G": {"-": 0.523076923, "T": 0.363636364, "C": 0.529411765, "A": 0.9},
            "U": {"-": 0.65625, "G": 0.625, "C": 0.647058824, "A": 0.476190476}},
        5: {"-": {"T": 0.133333333, "G": 0.071428571, "C": 0.24, "A": 0},
            "A": {"-": 0.328125, "T": 0.363636364, "G": 0.72, "C": 0.5},
            "C": {"-": 0.365079365, "T": 0.636363636, "G": 0.6, "A": 0.571428571},
            "G": {"-": 0.296875, "T": 0.3, "C": 0.785714286, "A": 0.866666667},
            "U": {"-": 0.215384615, "G": 0.64, "C": 1, "A": 0.5}},
        6: {"-": {"T": 0.066666667, "G": 0, "C": 0.142857143, "A": 0.083333333},
            "A": {"-": 0.276923077, "T": 0.714285714, "G": 0.714285714, "C": 0.454545455},
            "C": {"-": 0.4375, "T": 0.928571429, "G": 0.5, "A": 0.928571429},
            "G": {"-": 0.390625, "T": 0.666666667, "C": 0.681818182, "A": 1},
            "U": {"-": 0.446153846, "G": 0.571428571, "C": 0.909090909, "A": 0.866666667}},
        7: {"-": {"T": 0, "G": 0, "C": 0.058823529, "A": 0},
            "A": {"-": 0.430769231, "T": 0.4375, "G": 0.705882353, "C": 0.4375},
            "C": {"-": 0.569230769, "T": 0.8125, "G": 0.470588235, "A": 0.75},
            "G": {"-": 0.492063492, "T": 0.571428571, "C": 0.6875, "A": 1},
            "U": {"-": 0.553846154, "G": 0.588235294, "C": 0.6875, "A": 0.875}},
        8: {"-": {"T": 0, "G": 0, "C": 0.133333333, "A": 0.0625},
            "A": {"-": 0.234375, "T": 0.428571429, "G": 0.733333333, "C": 0.428571429},
            "C": {"-": 0.123076923, "T": 0.875, "G": 0.642857143, "A": 0.65},
            "G": {"-": 0.301587302, "T": 0.625, "C": 0.615384615, "A": 1},
            "U": {"-": 0.296875, "G": 0.733333333, "C": 1, "A": 0.8}},
        9: {"-": {"T": 0.142857143, "G": 0, "C": 0.238095238, "A": 0.1875},
            "A": {"-": 0.276923077, "T": 0.6, "G": 0.666666667, "C": 0.571428571},
            "C": {"-": 0.415384615, "T": 0.875, "G": 0.619047619, "A": 0.857142857},
            "G": {"-": 0.328125, "T": 0.533333333, "C": 0.538461538, "A": 0.642857143},
            "U": {"-": 0.609375, "G": 0.619047619, "C": 0.923076923, "A": 0.928571429}},
        10: {"-": {"T": 0.4, "G": 0.133333333, "C": 0.333333333, "A": 0.25},
             "A": {"-": 0.569230769, "T": 0.882352941, "G": 0.555555556, "C": 0.333333333},
             "C": {"-": 0.615384615, "T": 0.941176471, "G": 0.388888889, "A": 0.866666667},
             "G": {"-": 0.587301587, "T": 0.8125, "C": 0.4, "A": 0.933333333},
             "U": {"-": 0.650793651, "G": 0.5, "C": 0.533333333, "A": 0.857142857}},
        11: {"-": {"T": 0.3125, "G": 0.133333333, "C": 0.1, "A": 0.076923077},
             "A": {"-": 0.276923077, "T": 0.307692308, "G": 0.65, "C": 0.4},
             "C": {"-": 0.169230769, "T": 0.307692308, "G": 0.25, "A": 0.75},
             "G": {"-": 0.25, "T": 0.384615385, "C": 0.428571429, "A": 1},
             "U": {"-": 0.323076923, "G": 0.4, "C": 0.666666667, "A": 0.75}},
        12: {"-": {"T": 0.071428571, "G": 0, "C": 0.055555556, "A": 0.076923077},
             "A": {"-": 0.061538462, "T": 0.333333333, "G": 0.722222222, "C": 0.263157895},
             "C": {"-": 0.046875, "T": 0.538461538, "G": 0.444444444, "A": 0.714285714},
             "G": {"-": 0.125, "T": 0.384615385, "C": 0.529411765, "A": 0.933333333},
             "U": {"-": 0.030769231, "G": 0.5, "C": 0.947368421, "A": 0.8}},
        13: {"-": {"T": 0.076923077, "G": 0, "C": 0.043478261, "A": 0},
             "A": {"-": 0.03125, "T": 0.3, "G": 0.652173913, "C": 0.210526316},
             "C": {"-": 0, "T": 0.7, "G": 0.136363636, "A": 0.384615385},
             "G": {"-": 0.046875, "T": 0.3, "C": 0.421052632, "A": 0.923076923},
             "U": {"-": 0, "G": 0.260869565, "C": 0.789473684, "A": 0.692307692}},
        14: {"-": {"T": 0.047619048, "G": 0, "C": 0, "A": 0},
             "A": {"-": 0.046153846, "T": 0.533333333, "G": 0.466666667, "C": 0.214285714},
             "C": {"-": 0.046153846, "T": 0.733333333, "G": 0, "A": 0.35},
             "G": {"-": 0.063492063, "T": 0.266666667, "C": 0.428571429, "A": 0.75},
             "U": {"-": 0.061538462, "G": 0, "C": 0.285714286, "A": 0.619047619}},
        15: {"-": {"T": 0.052631579, "G": 0, "C": 0, "A": 0},
             "A": {"-": 0.03125, "T": 0.2, "G": 0.65, "C": 0.272727273},
             "C": {"-": 0, "T": 0.066666667, "G": 0.05, "A": 0.222222222},
             "G": {"-": 0, "T": 0.142857143, "C": 0.272727273, "A": 0.941176471},
             "U": {"-": 0, "G": 0.05, "C": 0.272727273, "A": 0.578947368}},
        16: {"-": {"T": 0, "G": 0, "C": 0, "A": 0},
             "A": {"-": 0.015384615, "T": 0, "G": 0.192307692, "C": 0},
             "C": {"-": 0, "T": 0.307692308, "G": 0.153846154, "A": 1},
             "G": {"-": 0, "T": 0, "C": 0, "A": 1},
             "U": {"-": 0, "G": 0.346153846, "C": 0.666666667, "A": 0.909090909}},
        17: {"-": {"T": 0, "G": 0.058823529, "C": 0, "A": 0},
             "A": {"-": 0, "T": 0.133333333, "G": 0.176470588, "C": 0.176470588},
             "C": {"-": 0, "T": 0.466666667, "G": 0.058823529, "A": 0.466666667},
             "G": {"-": 0, "T": 0.25, "C": 0.235294118, "A": 0.933333333},
             "U": {"-": 0, "G": 0.117647059, "C": 0.705882353, "A": 0.533333333}},
        18: {"-": {"T": 0.153846154, "G": 0.142857143, "C": 0.2, "A": 0},
             "A": {"-": 0, "T": 0.5, "G": 0.4, "C": 0.19047619},
             "C": {"-": 0, "T": 0.642857143, "G": 0.133333333, "A": 0.538461538},
             "G": {"-": 0, "T": 0.666666667, "C": 0.476190476, "A": 0.692307692},
             "U": {"-": 0, "G": 0.333333333, "C": 0.428571429, "A": 0.666666667}},
        19: {"-": {"T": 0.142857143, "G": 0.035714286, "C": 0.4375, "A": 0.307692308},
             "A": {"-": 0, "T": 0.538461538, "G": 0.375, "C": 0.206896552},
             "C": {"-": 0, "T": 0.461538462, "G": 0.125, "A": 0.428571429},
             "G": {"-": 0.063492063, "T": 0.666666667, "C": 0.448275862, "A": 0.714285714},
             "U": {"-": 0, "G": 0.25, "C": 0.275862069, "A": 0.285714286}},
        20: {"-": {"T": 0.6, "G": 0.045454545, "C": 0.529411765, "A": 0.6},
             "A": {"-": 0.061538462, "T": 0.6, "G": 0.764705882, "C": 0.227272727},
             "C": {"-": 0.03125, "T": 0.3, "G": 0.058823529, "A": 0.5},
             "G": {"-": 0.09375, "T": 0.7, "C": 0.428571429, "A": 0.9375},
             "U": {"-": 0, "G": 0.176470588, "C": 0.090909091, "A": 0.5625}}
    }

    @classmethod
    def off_target_score(cls, guide_sequence, target_sequence):
        """
        Calculates the score of a guide sequence and a target sequence using the class parameters for the
        endonuclease

        :param guide_sequence: str
            String corresponding to the guide RNA sequence (without the PAM)
        :param target_sequence: str
            String corresponding to the target site in the genome

        :return: score [numeric]
            Cutting frequency score of the target sequence with the specified guide sequence
        """

        guide_sequence, target_sequence, pam_seq = cls._process_sequences(guide_sequence, target_sequence)
        guide_sequence = str(guide_sequence).replace("T", "U")

        # Iterate through the PAM grid to find the cutting frequency of the PAM (default to 0)
        pam_array = cls.PAM_grid
        for i, p_chr in enumerate(list(pam_seq)):
            try:
                pam_array = pam_array[p_chr]
            except KeyError:
                return 0

        # Sanity check
        if isinstance(pam_array, numbers.Number):
            score = pam_array
        else:
            return 0

        # Walk through the guide sequence and compare it to the target sequence
        # If there is a mismatch, look up the mismatch score
        g_s = list(guide_sequence)
        t_s = list(target_sequence)

        for i, g_chr in enumerate(g_s):
            try:
                t_chr = t_s[i]
            except IndexError:
                t_chr = "-"

            if g_chr is t_chr:
                continue
            elif g_chr is "U" and t_chr is "T":
                continue

            try:
                score *= cls.spacer_grid[i + 1][g_chr][t_chr]
            except KeyError:
                return 0

        return score
