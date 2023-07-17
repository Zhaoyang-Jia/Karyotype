from Events import *

from itertools import groupby


def Output_KT(KT: {str: Chromosome}, output_path):
    # enumerate all segments in order
    def find_deleted_segments(arm: Arm):
        deleted_segments = []
        for event_itr in arm.history:
            if event_itr[0] == 'deletion':
                deleted_segments.extend(event_itr[1])
        return deleted_segments

    def update_set_with_positive_direction(working_set, segments_to_add):
        for segment_itr in segments_to_add:
            if segment_itr.direction():
                working_set.add(segment_itr)
            else:
                new_segment = segment_itr.duplicate()
                new_segment.invert()
                working_set.add(new_segment)
        return working_set

    segment_set = set()
    for chr_name, chr_obj in KT.items():
        segment_set = update_set_with_positive_direction(segment_set, chr_obj.p_arm.segments)
        segment_set = update_set_with_positive_direction(segment_set, chr_obj.q_arm.segments)
        segment_set = update_set_with_positive_direction(segment_set, find_deleted_segments(chr_obj.p_arm))
        segment_set = update_set_with_positive_direction(segment_set, find_deleted_segments(chr_obj.q_arm))
    segment_list = list(segment_set)
    segment_list.sort()

    # segment_list = []
    # for chr_name, chr_obj in KT.items():
    #     segment_list.extend(chr_obj.p_arm.segments)
    #     segment_list.extend(chr_obj.q_arm.segments)
    #     segment_list.extend(find_deleted_segments(chr_obj.p_arm))
    #     segment_list.extend(find_deleted_segments(chr_obj.q_arm))
    # segment_list.sort()

    # put segments into right direction
    # for segment in segment_list:
    #     if not segment.direction():
    #         segment.invert()
    #
    # key_function = lambda x: (x.chr, x.start, x.end)
    # unique_segments = [next(group) for key, group in groupby(segment_list, key_function)]
    # segment_list = unique_segments

    segment_dict = {}
    for index in range(len(segment_list)):
        segment_dict[segment_list[index]] = index + 1

    # enumerate event history
    def output_history(arm: Arm):
        history_list = []
        for history_itr in arm.history:
            full_event_name = history_itr[0]
            event_segments = history_itr[1]

            indexed_segments = []
            for segment_itr in event_segments:
                for segment_dict_segment in segment_dict:
                    if segment_itr.same_segment_ignore_dir(segment_dict_segment):
                        indexed_segments.append(str(segment_dict[segment_dict_segment]))
                        break
            history_list.append([full_event_name, indexed_segments])
        return history_list

    event_history = {}  # {chr_name: [[event_name, [segment_index]]]}
    for chr_name, chr_obj in KT.items():
        current_chr_history = output_history(chr_obj.p_arm)
        current_chr_history.extend(output_history(chr_obj.q_arm))
        if len(current_chr_history) == 0:
            continue
        else:
            event_history[chr_name] = current_chr_history

    # enumerate chromosome KT
    def naming_segments(arm: Arm):
        str_list = []
        for segment_itr in arm.segments:
            if segment_itr.direction():
                segment_direction = "+"
                positive_direction_representation = segment_itr
            else:
                segment_direction = "-"
                new_segment = segment_itr.duplicate()
                new_segment.invert()
                positive_direction_representation = new_segment
            segment_sorted_index = segment_dict[positive_direction_representation]
            str_list.append(str(segment_sorted_index) + segment_direction)
        return str_list

    chromosome_KT = {}  # {chr_name: ([segment_index], [segment_index])}
    for chr_name, chr_obj in KT.items():
        p_str = naming_segments(chr_obj.p_arm)
        q_str = naming_segments(chr_obj.q_arm)

        chromosome_KT[chr_name] = (p_str, q_str)

    # output
    with open(output_path, 'w') as fp_write:
        # enumerate segments
        for segment, index in segment_dict.items():
            fp_write.writelines("{index}\t{segment_str}\n"
                                .format(index=index, segment_str=segment.to_string_ignore_dir()))

        # write history
        fp_write.writelines("---\n")
        for chr_name, chr_history in event_history.items():
            fp_write.writelines("{}:\n".format(chr_name))
            for event in chr_history:
                fp_write.writelines("{} on segment {}\n".format(event[0], ','.join(event[1])))

        # write KT
        fp_write.writelines("---\n")
        for chromosome_name, chromosome_segments in chromosome_KT.items():
            fp_write.writelines("{chr_name}: {p_str}\t{q_str}\n"
                                .format(chr_name=chromosome_name,
                                        p_str=','.join(chromosome_segments[0]),
                                        q_str=','.join(chromosome_segments[1])))


test_arm1 = Arm([Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
                 Segment('Chr1', 38, 39), Segment('Chr1', 40, 50), Segment('Chr1', 51, 76)])
test_arm2 = Arm([Segment('Chr1', 100, 200)])
test_centromere = Arm([Segment('Chr1', 77, 99)])
Chr1 = Chromosome('Chr1', test_arm1, test_arm2, 100, 100, test_centromere)

test_arm3 = Arm([Segment('Chr2', 0, 25), Segment('Chr2', 26, 37),
                 Segment('Chr2', 38, 39), Segment('Chr2', 40, 50), Segment('Chr2', 51, 76)])
test_arm4 = Arm([Segment('Chr2', 100, 200)])
test_centromere = Arm([Segment('Chr2', 77, 99)])
Chr2 = Chromosome('Chr2', test_arm3, test_arm4, 100, 100, test_centromere)

this_KT = {'Chr1': Chr1, 'Chr2': Chr2}
# this_KT = {'Chr1': Chr1}
# duplication(this_KT['Chr1'].p_arm, 27, 53)
# duplication_inversion(this_KT['Chr1'].p_arm, 27, 53)
# translocation_reciprocal(this_KT['Chr1'].p_arm, 27, 53, this_KT['Chr2'].q_arm, 27, 53)
# inversion(this_KT['Chr1'].p_arm, 27, 53)
# deletion(this_KT['Chr2'].q_arm, 27, 53)

Output_KT(this_KT, "test_KT_out4.txt")
