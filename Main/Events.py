from Prepare_Raw_KT import Prepare_Raw_KT


def Segment_Length(start, end):
    return abs(start - end + 1)


def Segment_Direction(start, end):
    """
    :return: 1 for +, 0 for -
    """
    return start <= end


def Segment_Indexing(KT_arm):
    segment_indices = []
    current_index = 0
    for segment in KT_arm:
        next_index = current_index + Segment_Length(segment[1], segment[2]) - 1
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
                segment_indices[segment_index][0] = cut_low + 1
                break
            elif cut_high >= segment_high:
                # complete deletion, continue
                index_marked_for_removal.append(segment_index)
                cut_low = segment_high + 1
        elif cut_low > segment_low:
            if cut_high >= segment_high:
                # right deletion, continue
                segment_indices[segment_index][1] = cut_low - 1
                cut_low = segment_high + 1
            elif cut_high < segment_high:
                # middle deletion, split segment, terminate
                segment_indices[segment_index][1] = cut_low - 1
                new_segment = (current_arm[segment_index])


KT = Prepare_Raw_KT(['Chr1', 'ChrX'], "../Metadata/Full_Genome_Indices.txt")
print(Segment_Indexing(KT[0]['p'])[0][1])
