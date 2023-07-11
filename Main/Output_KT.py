from Events import *


def Output_KT(KT, output_path):
    # enumerate all segments while outputting KT
    segment_list = []
    enumerate_index = 1
    KT_lines = []
    for chromosome in KT:
        p_line = []
        q_line = []
        for p_segment in chromosome['p']:
            segment_list.append(p_segment)
            if p_segment.direction():
                p_line.append("+" + str(enumerate_index))
            else:
                p_line.append("-" + str(enumerate_index))
            enumerate_index += 1
        for q_segment in chromosome['q']:
            segment_list.append(q_segment)
            if q_segment.direction():
                q_line.append("+" + str(enumerate_index))
            else:
                q_line.append("-" + str(enumerate_index))
            enumerate_index += 1

        new_KT_line = ",".join(p_line) + "\t" + ",".join(q_line)
        KT_lines.append(new_KT_line)

    # output
    with open(output_path, 'w') as fp_write:
        # enumerate segments
        for index in range(len(segment_list)):
            current_segment = segment_list[index]
            ordered_start = min(current_segment.start, current_segment.end)
            ordered_end = max(current_segment.start, current_segment.end)
            fp_write.writelines("{enumerator}\t{chr}:{start}-{end}\n"
                                .format(enumerator=index + 1,
                                        chr=current_segment.chr,
                                        start=ordered_start,
                                        end=ordered_end))
        # write KT
        fp_write.writelines("\n")
        fp_write.writelines("\n".join(KT_lines))


test_p_arm1 = [Segment('Chr1', 0, 25), Segment('Chr1', 26, 37),
               Segment('Chr1', 38, 39), Segment('Chr1', 40, 50), Segment('Chr1', 51, 76)]
test_p_arm2 = [Segment('Chr1', 0, 100)]
# this_KT = Prepare_Raw_KT(['Chr1'], "../Metadata/Full_Genome_Indices.txt")
this_KT = [{'p': test_p_arm1, 'q': []}]
inversion(this_KT, 0, 'p', 27, 53)

Output_KT(this_KT, "test_KT_out.txt")