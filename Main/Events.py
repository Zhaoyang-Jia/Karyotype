class Segment:
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end

    def __len__(self):
        return abs(self.end - self.start) + 1

    def __str__(self):
        return "({}, {}, {})".format(self.chr, self.start, self.end)

    def direction(self):
        """
        :return: 1 for +, 0 for -
        """
        return self.start <= self.end

    def duplicate(self):
        return Segment(self.chr, self.start, self.end)

    def left_deletion(self, end_index):
        """
        :param end_index: grounded index, right-bounded
        :return: None
        """
        if self.direction():
            self.start = self.start + end_index + 1
        else:
            self.start = self.start - end_index - 1

    def right_deletion(self, start_index):
        """
        :param start_index: grounded index, left-bounded
        :return: None
        """
        if self.direction():
            self.end = self.start + start_index - 1
        else:
            self.end = self.start - start_index + 1

    def inversion(self):
        temp_start = self.start
        self.start = self.end
        self.end = temp_start


def Prepare_Raw_KT(chr_of_interest, genome_index_file):
    """
    Compose unedited KT
    :param chr_of_interest: list of chromosomes to generate KT
    :param genome_index_file: .txt metadata file containing centromere, telomere, and genome length information
    :return:
    """
    chromosomes = []
    with open(genome_index_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            if line[0] in chr_of_interest:
                new_chromosome = {"p": [], "q": []}
                chromosome_origin = line[0]

                # p-arm
                start_index = int(line[2])
                end_index = int(line[3])
                new_chromosome['p'].append(Segment(chromosome_origin, start_index, end_index))

                # q-arm
                start_index = int(line[4])
                end_index = int(line[5])
                new_chromosome['q'].append(Segment(chromosome_origin, start_index, end_index))
                chromosomes.append(new_chromosome)

    return chromosomes


def Segment_Indexing(KT_arm):
    segment_indices = []
    current_index = 0
    for segment in KT_arm:
        next_index = current_index + len(segment) - 1
        segment_indices.append([current_index, next_index])
        current_index = next_index + 1
    return segment_indices


def deletion(KT, chromosome_index, arm, cut_low, cut_high):
    current_arm = KT[chromosome_index][arm]
    segment_indices = Segment_Indexing(current_arm)

    index_marked_for_removal = []
    for segment_index in range(len(current_arm)):
        segment_low = segment_indices[segment_index][0]
        segment_high = segment_indices[segment_index][1]
        if cut_low == segment_low:
            if cut_high < segment_high:
                # left deletion, terminate
                grounded_index = cut_high - segment_low
                current_arm[segment_index].left_deletion(grounded_index)
                break
            elif cut_high >= segment_high:
                # complete deletion, continue
                index_marked_for_removal.append(segment_index)
                cut_low = segment_high + 1
        elif segment_low < cut_low <= segment_high:
            if cut_high >= segment_high:
                # right deletion, continue
                grounded_index = cut_low - segment_low
                current_arm[segment_index].right_deletion(grounded_index)
                cut_low = segment_high + 1
            elif cut_high < segment_high:
                # middle deletion, split segment, terminate
                new_segment = current_arm[segment_index].duplicate()
                left_grounded_index = cut_low - segment_low
                right_grounded_index = cut_high - segment_low
                current_arm[segment_index].right_deletion(left_grounded_index)
                new_segment.left_deletion(right_grounded_index)
                current_arm.insert(segment_index + 1, new_segment)
                break

    # remove empty segments
    KT[chromosome_index][arm] = [element for index, element in enumerate(current_arm) if index not in index_marked_for_removal]


test_p_arm1 = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 38), Segment('Chr1', 39, 50), Segment('Chr1', 51, 76)]
test_p_arm2 = [Segment('Chr1', 0, 100)]
# this_KT = Prepare_Raw_KT(['Chr1'], "../Metadata/Full_Genome_Indices.txt")
this_KT = [{'p': test_p_arm2}]
deletion(this_KT, 0, 'p', 27, 53)

print(this_KT[0]['p'][0])
