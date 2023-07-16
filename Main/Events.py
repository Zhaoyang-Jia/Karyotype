import copy


class Segment:
    def __init__(self, chr: str, start: int, end: int):
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

    def left_delete(self, bp_to_delete):
        """
        :param bp_to_delete: number of bp deleting
        :return: None
        """
        if self.direction():
            self.start = self.start + bp_to_delete
        else:
            self.start = self.start - bp_to_delete

    def right_delete(self, bp_to_delete):
        """
        :param bp_to_delete: number of bp deleting
        :return: None
        """
        if self.direction():
            self.end = self.end - bp_to_delete
        else:
            self.end = self.start + bp_to_delete

    def invert(self):
        temp_start = self.start
        self.start = self.end
        self.end = temp_start


class Arm:
    def __init__(self, segments: [Segment]):
        self.segments = segments

    def __len__(self):
        current_sum = 0
        for segment in self.segments:
            current_sum += len(segment)
        return current_sum

    def __str__(self):
        return_str = ''
        for segment in self.segments:
            return_str += str(segment) + '\n'
        return return_str

    def generate_breakpoint(self, breakpoint_index: int):
        """
        split segment such that the breakpoint_index is garenteed to be the end index of a Segment
        :param breakpoint_index:
        :return: None
        """
        current_segment_index = 0
        current_bp_index = -1  # corrects 0-index off-shift
        for segment in self.segments:
            current_bp_index += len(segment)
            if current_bp_index == breakpoint_index:
                # breakpoint exists
                return
            elif current_bp_index > breakpoint_index:
                # breakpoint within the current segment
                previous_bp_index = current_bp_index - len(segment)
                new_segment = segment.duplicate()
                new_segment.left_delete(breakpoint_index - previous_bp_index)
                segment.right_delete(current_bp_index - breakpoint_index)
                self.segments.insert(current_segment_index + 1, new_segment)
                return
            else:
                # breakpoint location not yet met
                current_segment_index += 1

    def delete_segment_by_indices(self, indices):
        self.segments = [element for index, element in enumerate(self.segments) if index not in indices]


class Chromosome:
    def __init__(self, name: str, p_arm: Arm, q_arm: Arm, t1_len: int, t2_len: int, centromere: Arm):
        self.name = name
        self.p_arm = p_arm
        self.q_arm = q_arm
        self.t1_len = t1_len
        self.t2_len = t2_len
        self.centromere = centromere

    def p_arm_len(self):
        return len(self.p_arm)

    def q_arm_len(self):
        return len(self.q_arm)


def Prepare_Raw_KT(chr_of_interest: [str], genome_index_file: str) -> [Chromosome]:
    """
    Compose unedited KT
    :param chr_of_interest: list of chromosomes to generate KT
    :param genome_index_file: .txt metadata file containing centromere, telomere, and genome length information
    :return: a list of Chromosome Object
    """
    chromosomes = []
    with open(genome_index_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            if line[0] in chr_of_interest:
                chr_name = line[0]
                p_arm_segment = Segment(chr_name, int(line[2]), int(line[3]))
                q_arm_segment = Segment(chr_name, int(line[4]), int(line[5]))
                t1_len = int(line[2])
                t2_len = int(line[1]) - int(line[5]) - 1
                centromere_segment = Segment(chr_name, int(line[3]) + 1, int(line[4]) - 1)

                chromosomes.append(Chromosome(chr_name, Arm([p_arm_segment]), Arm([q_arm_segment]),
                                              t1_len, t2_len, Arm([centromere_segment])))

    return chromosomes


# def Segment_Indexing(KT_arm):
#     segment_indices = []
#     current_index = 0
#     for segment in KT_arm:
#         next_index = current_index + len(segment) - 1
#         segment_indices.append([current_index, next_index])
#         current_index = next_index + 1
#     return segment_indices


def deletion(event_arm: Arm, left_event_index: int, right_event_index: int) -> [Segment]:
    """
    perform deletion event, inplace
    :param event_arm: chromosome arm that the event will happen in
    :param left_event_index: beginning of deletion, this index will be deleted
    :param right_event_index: end of deletion, this index will be deleted
    :return: list of Segments that got deleted
    """
    event_arm.generate_breakpoint(left_event_index - 1)
    event_arm.generate_breakpoint(right_event_index)

    segment_index_for_deletion = []
    current_segment_index = 0
    current_bp_index = -1  # corrects 0-index off-shift
    for segment in event_arm.segments:
        current_bp_index += len(segment)
        if left_event_index <= current_bp_index <= right_event_index:
            segment_index_for_deletion.append(current_segment_index)
        elif current_bp_index > right_event_index:
            break
        current_segment_index += 1

    # document segments deleted
    segments_deleted = []
    for segment_index in range(len(event_arm.segments)):
        if segment_index in segment_index_for_deletion:
            document_segment = event_arm.segments[segment_index].duplicate()
            segments_deleted.append(document_segment)
    # remove empty segments
    event_arm.delete_segment_by_indices(segment_index_for_deletion)

    return segments_deleted


# def duplication(KT, chromosome_index, arm, cut_low, cut_high):
#     current_arm = KT[chromosome_index][arm]
#
#     segment_indices = Segment_Indexing(current_arm)
#     # split left boundary into two segments
#     for segment_index in range(len(current_arm)):
#         segment_low = segment_indices[segment_index][0]
#         segment_high = segment_indices[segment_index][1]
#         if cut_low == segment_low:
#             # then no split required
#             # note, if cut_low == segment_high, split is still required
#             break
#         if segment_low < cut_low <= segment_high:
#             # split
#             new_segment = current_arm[segment_index].duplicate()
#             right_deletion_grounded_index = cut_low - segment_low
#             left_deletion_grounded_index = cut_low - segment_low - 1
#             current_arm[segment_index].right_deletion(right_deletion_grounded_index)
#             new_segment.left_deletion(left_deletion_grounded_index)
#             current_arm.insert(segment_index + 1, new_segment)
#             break
#
#     segment_indices = Segment_Indexing(current_arm)
#     # split right boundary into two segments
#     for segment_index in range(len(current_arm)):
#         segment_low = segment_indices[segment_index][0]
#         segment_high = segment_indices[segment_index][1]
#         if cut_high == segment_high:
#             # then no split required
#             # note, if cut_high == segment_low, split is still required
#             break
#         if segment_low <= cut_high < segment_high:
#             # split
#             new_segment = current_arm[segment_index].duplicate()
#             right_deletion_grounded_index = cut_high - segment_low + 1
#             left_deletion_grounded_index = cut_high - segment_low
#             current_arm[segment_index].right_deletion(right_deletion_grounded_index)
#             new_segment.left_deletion(left_deletion_grounded_index)
#             current_arm.insert(segment_index + 1, new_segment)
#             break
#
#     segment_indices = Segment_Indexing(current_arm)
#     # duplicate segments
#     duplicated_segments = []
#     first_index = -1
#     for segment_index in range(len(current_arm)):
#         segment_low = segment_indices[segment_index][0]
#         segment_high = segment_indices[segment_index][1]
#         if cut_low <= segment_low <= segment_high <= cut_high:
#             if first_index == -1:
#                 first_index = segment_index
#             duplicated_segments.append(current_arm[segment_index].duplicate())
#
#     # duplicate right before the first segment for duplication
#     KT[chromosome_index][arm][first_index :first_index] = duplicated_segments
#
#
# def inversion(KT, chromosome_index, arm, cut_low, cut_high):
#     current_arm = KT[chromosome_index][arm]
#
#     segment_indices = Segment_Indexing(current_arm)
#     # split left boundary into two segments
#     for segment_index in range(len(current_arm)):
#         segment_low = segment_indices[segment_index][0]
#         segment_high = segment_indices[segment_index][1]
#         if cut_low == segment_low:
#             # then no split required
#             # note, if cut_low == segment_high, split is still required
#             break
#         if segment_low < cut_low <= segment_high:
#             # split
#             new_segment = current_arm[segment_index].duplicate()
#             right_deletion_grounded_index = cut_low - segment_low
#             left_deletion_grounded_index = cut_low - segment_low - 1
#             current_arm[segment_index].right_deletion(right_deletion_grounded_index)
#             new_segment.left_deletion(left_deletion_grounded_index)
#             current_arm.insert(segment_index + 1, new_segment)
#             break
#
#     segment_indices = Segment_Indexing(current_arm)
#     # split right boundary into two segments
#     for segment_index in range(len(current_arm)):
#         segment_low = segment_indices[segment_index][0]
#         segment_high = segment_indices[segment_index][1]
#         if cut_high == segment_high:
#             # then no split required
#             # note, if cut_high == segment_low, split is still required
#             break
#         if segment_low <= cut_high < segment_high:
#             # split
#             new_segment = current_arm[segment_index].duplicate()
#             right_deletion_grounded_index = cut_high - segment_low + 1
#             left_deletion_grounded_index = cut_high - segment_low
#             current_arm[segment_index].right_deletion(right_deletion_grounded_index)
#             new_segment.left_deletion(left_deletion_grounded_index)
#             current_arm.insert(segment_index + 1, new_segment)
#             break
#
#     segment_indices = Segment_Indexing(current_arm)
#     # mark segments for inversion
#     segment_index_to_invert = []
#     for segment_index in range(len(current_arm)):
#         segment_low = segment_indices[segment_index][0]
#         segment_high = segment_indices[segment_index][1]
#         if cut_low <= segment_low <= segment_high <= cut_high:
#             segment_index_to_invert.append(segment_index)
#
#     # invert
#     new_arm = []
#     first_index_for_inversion = segment_index_to_invert[0]
#     for segment_index in range(len(current_arm)):
#         if segment_index < first_index_for_inversion:
#             new_arm.append(current_arm[segment_index])
#         else:
#             break
#     for segment_index in reversed(segment_index_to_invert):
#         current_arm[segment_index].invert()
#         new_arm.append(current_arm[segment_index])
#     for segment_index in range(segment_index_to_invert[-1] + 1, len(current_arm)):
#         new_arm.append(current_arm[segment_index])
#
#     KT[chromosome_index][arm] = new_arm


# test_p_arm1 = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
#                Segment('Chr1', 38, 39), Segment('Chr1', 40, 50), Segment('Chr1', 51, 76)]
# test_p_arm2 = [Segment('Chr1', 0, 100)]
# # this_KT = Prepare_Raw_KT(['Chr1'], "../Metadata/Full_Genome_Indices.txt")
# this_Arm = Arm(test_p_arm1)
# documentation = deletion(this_Arm, 27, 53)
# print(this_Arm)
# print(Arm(documentation))
