import copy

from Structures import *


def Prepare_Raw_KT(chr_of_interest: [str], genome_index_file: str) -> {str: Chromosome}:
    """
    Compose unedited KT
    :param chr_of_interest: list of chromosomes to generate KT
    :param genome_index_file: .txt metadata file containing centromere, telomere, and genome length information
    :return: a list of Chromosome Object
    """
    chromosomes = {}
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

                chromosomes[chr_name] = Chromosome(chr_name, Arm([p_arm_segment]), Arm([q_arm_segment]),
                                                   t1_len, t2_len, Arm([centromere_segment]))

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
    segments_selected_indices = []
    current_segment_index = 0
    current_bp_index = -1  # corrects 0-index off-shift
    for segment in event_arm.segments:
        current_bp_index += len(segment)
        if left_event_index <= current_bp_index <= right_event_index:
            segments_selected.append(segment)
            segments_selected_indices.append(current_segment_index)
        elif current_bp_index > right_event_index:
            break
        current_segment_index += 1

    return segments_selected, segments_selected_indices


def deletion(event_arm: Arm, left_event_index: int, right_event_index: int):
    """
    perform deletion event, inplace
    :param event_arm: chromosome arm that the event will happen in
    :param left_event_index: beginning of deletion, this index will be deleted
    :param right_event_index: end of deletion, this index will be deleted
    :return: None
    """
    event_segments, event_segments_indices = locate_segments_for_event(event_arm, left_event_index, right_event_index)
    # document segments deleted
    event_arm.append_history('deletion', event_segments)
    # remove empty segments
    event_arm.delete_segments_by_index(event_segments_indices)

    return


def duplication(event_arm: Arm, left_event_index: int, right_event_index: int):
    """
    duplication even, inplace
    :param event_arm: chromosome arm that the event will happen in
    :param left_event_index: beginning of deletion, this index will be deleted
    :param right_event_index: end of deletion, this index will be deleted
    :return: None
    """
    event_segments, event_segment_indices = locate_segments_for_event(event_arm, left_event_index, right_event_index)
    # document segments duplicated
    event_arm.append_history('duplication', event_segments)
    # duplicate segments
    event_arm.duplicate_segments_by_index(event_segment_indices)


def inversion(event_arm: Arm, left_event_index: int, right_event_index: int):
    """
    inversion even, inplace
    :param event_arm: chromosome arm tha t the event will happen in
    :param left_event_index: beginning of deletion, this index will be deleted
    :param right_event_index: end of deletion, this index will be deleted
    :return: None
    """
    event_segments, event_segments_indices = locate_segments_for_event(event_arm, left_event_index, right_event_index)
    # document segments inverted
    event_arm.append_history('inversion', event_segments)
    # invert segments
    event_arm.invert_segments_by_index(event_segments_indices)


def translocation_reciprocal(event_arm1: Arm, arm1_left_index: int, arm1_right_index: int,
                             event_arm2: Arm, arm2_left_index: int, arm2_right_index: int):
    arm1_segments, arm1_segment_indices = locate_segments_for_event(event_arm1, arm1_left_index, arm1_right_index)
    arm2_segments, arm2_segment_indices = locate_segments_for_event(event_arm2, arm2_left_index, arm2_right_index)
    arm1_start_segment_index = arm1_segment_indices[0]
    arm2_start_segment_index = arm2_segment_indices[0]

    event_arm1.append_history('reciprocal translocation', arm1_segments)
    event_arm2.append_history('reciprocal translocation', arm2_segments)

    event_arm1.delete_segments_by_index(arm1_segment_indices)
    event_arm2.delete_segments_by_index(arm2_segment_indices)
    event_arm2.segments[arm2_start_segment_index:arm2_start_segment_index] = arm1_segments
    event_arm1.segments[arm1_start_segment_index:arm1_start_segment_index] = arm2_segments


def duplication_inversion(event_arm: Arm, left_event_index: int, right_event_index: int):
    event_segments, event_segment_indices = locate_segments_for_event(event_arm, left_event_index, right_event_index)
    event_arm.append_history('duplication inversion', event_segments)
    new_segment_start_index = event_segment_indices[-1] + 1
    new_segment_end_index = new_segment_start_index + len(event_segments) - 1
    event_arm.duplicate_segments_by_index(event_segment_indices)
    segments_for_inversion_indices = range(new_segment_start_index, new_segment_end_index + 1)
    event_arm.invert_segments_by_index(segments_for_inversion_indices)


# test_p_arm1 = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
#                Segment('Chr1', 38, 39), Segment('Chr1', 40, 50), Segment('Chr1', 51, 76)]
# test_p_arm2 = [Segment('Chr1', 0, 100)]
# # this_KT = Prepare_Raw_KT(['Chr1'], "../Metadata/Full_Genome_Indices.txt")
# Arm1 = Arm(test_p_arm1)
# Arm2 = Arm(test_p_arm2)
# duplication(Arm2, 27, 53)
# deletion(Arm2, 27, 51)
# print(Arm2)
