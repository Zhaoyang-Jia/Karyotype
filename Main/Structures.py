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

    def __lt__(self, other):
        def get_chr_order(chromosome_name):
            chr_extracted = chromosome_name.replace('Chr', '')
            if chr_extracted == 'X':
                return 23
            elif chr_extracted == 'Y':
                return 24
            else:
                return int(chr_extracted)

        if get_chr_order(self.chr) < get_chr_order(other.chr):
            return True
        elif get_chr_order(self.chr) > get_chr_order(other.chr):
            return False
        elif get_chr_order(self.chr) == get_chr_order(other.chr):
            return max(self.start, self.end) < max(other.start, other.end)

    def __eq__(self, other):
        if isinstance(other, Segment):
            return (self.chr, self.start, self.end) == (other.chr, other.start, other.end)
        return False

    def __hash__(self):
        return hash((self.chr, self.start, self.end))

    def __str__(self):
        return "({}, {}, {})".format(self.chr, self.start, self.end)

    def same_segment_ignore_dir(self, other):
        if self.start != other.start and self.start != other.end:
            return False
        if self.end != other.start and self.end != other.end:
            return False
        return True

    def to_string_ignore_dir(self):
        if self.direction():
            return "{}:{}-{}".format(self.chr, self.start, self.end)
        else:
            return "{}:{}-{}".format(self.chr, self.end, self.start)

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

    def delete_segments_by_index(self, segment_indices):
        self.segments = [segment for index, segment in enumerate(self.segments) if index not in segment_indices]

    def duplicate_segments_by_index(self, segment_indices):
        # segments come in order, so insert before the first segment
        index_of_insertion = segment_indices[0]
        new_segments = []
        for index in segment_indices:
            new_segments.append(self.segments[index].duplicate())

        self.segments[index_of_insertion:index_of_insertion] = new_segments

    def invert_segments_by_index(self, segment_indices):
        # segments come in order
        index_of_insertion = segment_indices[0]
        new_segments = []
        for index in reversed(segment_indices):
            new_segment = self.segments[index].duplicate()
            new_segment.invert()
            new_segments.append(new_segment)
        self.delete_segments_by_index(segment_indices)
        self.segments[index_of_insertion:index_of_insertion] = new_segments


class Chromosome:
    name: str
    p_arm: Arm
    q_arm: Arm
    centromere: Arm
    t1_len: int
    t2_len: int

    def __init__(self, name: str, p_arm: Arm, q_arm: Arm, t1_len: int, t2_len: int, centromere: Arm):
        self.name = name
        self.p_arm = p_arm
        self.q_arm = q_arm
        self.t1_len = t1_len
        self.t2_len = t2_len
        self.centromere = centromere

    def __len__(self):
        return self.p_arm_len() + self.q_arm_len()

    def p_arm_len(self):
        return len(self.p_arm)

    def q_arm_len(self):
        return len(self.q_arm)


class Genome:
    full_KT: {str: [Chromosome]}  # has exactly 24 slots, corresponding to the 24 possible chromosome types
    history: [(str, [Segment], Chromosome, Chromosome)]  # event type, event segments, chr from, chr to

    def __init__(self):
        # construct KT slots
        full_KT_keys = [str(i) for i in range(23)]
        full_KT_keys.extend(['X', 'Y'])
        for key in full_KT_keys:
            chr_name = 'Chr' + key
            self.full_KT[chr_name] = []
