import copy


class Segment:
    chr: str
    start: int
    end: int

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
    segments: [Segment]
    history: [(str, [Segment])]

    def __init__(self, segments: [Segment]):
        self.segments = segments
        self.history = []

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

    def append_history(self, event_type: str, segments: [Segment]):
        """
        add the newest event to history log
        :param event_type: del, inv, dup
        :param segments: segments that were selected for the event (e.g. if +8 inverted to -8, +8 is recorded)
        :return: None
        """
        new_history = tuple([event_type, segments])
        self.history.append(new_history)

    def delete_segments(self, segments_to_delete):
        self.segments = [obj for obj in self.segments if obj not in segments_to_delete]

    def duplicate_segments(self, segments_to_duplicate):
        # segments come in order, so insert before the first segment
        index_of_insertion = self.segments.index(segments_to_duplicate[0])
        new_segments = []
        for segment in segments_to_duplicate:
            new_segments.append(segment.duplicate())

        self.segments[index_of_insertion:index_of_insertion] = new_segments

    def invert_segments(self, segment_to_invert):
        # segments come in order
        index_of_insertion = self.segments.index(segment_to_invert[0])
        new_segments = []
        for segment in reversed(segment_to_invert):
            new_segment = segment.duplicate()
            new_segment.invert()
            new_segments.append(new_segment)
        self.delete_segments(segment_to_invert)
        self.segments[index_of_insertion:index_of_insertion] = new_segments


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


def locate_segments_for_event(event_arm: Arm, left_event_index: int, right_event_index: int) -> [Segment]:
    """
    create breakpoint and select the Segement between the breakpoints
    :param event_arm: chromosome arm that the event will happen in
    :param left_event_index: beginning of deletion, this index will be deleted
    :param right_event_index: end of deletion, this index will be deleted
    :return: a list of Segment for processing the event
    """
    event_arm.generate_breakpoint(left_event_index - 1)
    event_arm.generate_breakpoint(right_event_index)

    segments_selected = []
    current_bp_index = -1  # corrects 0-index off-shift
    for segment in event_arm.segments:
        current_bp_index += len(segment)
        if left_event_index <= current_bp_index <= right_event_index:
            segments_selected.append(segment)
        elif current_bp_index > right_event_index:
            break

    return segments_selected


def deletion(event_arm: Arm, left_event_index: int, right_event_index: int):
    """
    perform deletion event, inplace
    :param event_arm: chromosome arm that the event will happen in
    :param left_event_index: beginning of deletion, this index will be deleted
    :param right_event_index: end of deletion, this index will be deleted
    :return: None
    """
    segments_for_deletion = locate_segments_for_event(event_arm, left_event_index, right_event_index)
    # document segments deleted
    event_arm.append_history('del', segments_for_deletion)
    # remove empty segments
    event_arm.delete_segments(segments_for_deletion)

    return


def duplication(event_arm: Arm, left_event_index: int, right_event_index: int):
    """
    duplication even, inplace
    :param event_arm: chromosome arm that the event will happen in
    :param left_event_index: beginning of deletion, this index will be deleted
    :param right_event_index: end of deletion, this index will be deleted
    :return: None
    """
    segments_for_duplication = locate_segments_for_event(event_arm, left_event_index, right_event_index)
    # document segments duplicated
    event_arm.append_history('dup', segments_for_duplication)
    # duplicate segments
    event_arm.duplicate_segments(segments_for_duplication)


def inversion(event_arm: Arm, left_event_index: int, right_event_index: int):
    segments_for_inversion = locate_segments_for_event(event_arm, left_event_index, right_event_index)
    # document segments inverted
    event_arm.append_history('inv', segments_for_inversion)
    # invert segments
    event_arm.invert_segments(segments_for_inversion)


# test_p_arm1 = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
#                Segment('Chr1', 38, 39), Segment('Chr1', 40, 50), Segment('Chr1', 51, 76)]
# test_p_arm2 = [Segment('Chr1', 0, 100)]
# # this_KT = Prepare_Raw_KT(['Chr1'], "../Metadata/Full_Genome_Indices.txt")
# Arm1 = Arm(test_p_arm1)
# Arm2 = Arm(test_p_arm2)
# inversion(Arm1, 27, 53)
# print(Arm1)

