from Events import *


def Output_KT(KT: {str: Chromosome}, output_path):
    # enumerate all segments in order
    def find_deleted_segments(arm: Arm):
        deleted_segments = []
        for event in arm.history:
            if event[0] == 'del':
                deleted_segments.extend(event[1])
        return deleted_segments

    segment_list = []
    for chr_name, chr_obj in KT.items():
        segment_list.extend(chr_obj.p_arm.segments)
        segment_list.extend(chr_obj.q_arm.segments)
        segment_list.extend(find_deleted_segments(chr_obj.p_arm))
        segment_list.extend(find_deleted_segments(chr_obj.q_arm))
    segment_list.sort()

    segment_dict = {}
    for index in range(len(segment_list)):
        segment_dict[segment_list[index]] = index + 1

    # enumerate event history
    def output_history(arm: Arm):
        history_list = []
        for history_itr in arm.history:
            event_name = history_itr[0]
            event_segments = history_itr[1]
            full_event_name = ''
            if event_name == 'del':
                full_event_name = 'deletion'
            elif event_name == 'inv':
                full_event_name = 'inversion'
            elif event_name == 'dup:':
                full_event_name = 'dup'
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
            segment_sorted_index = segment_dict[segment_itr]
            if segment_itr.direction():
                segment_direction = "+"
            else:
                segment_direction = "-"
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

# this_KT = Prepare_Raw_KT(['Chr1'], "../Metadata/Full_Genome_Indices.txt")
this_KT = {'Chr1': Chr1, 'Chr2': Chr2}
inversion(this_KT['Chr1'].p_arm, 27, 53)
deletion(this_KT['Chr2'].q_arm, 27, 53)

Output_KT(this_KT, "test_KT_out2.txt")
