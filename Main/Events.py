def Parse_KT(KT_file):
    """
    Parse KT into p-arms and q-arms
    :param KT_file:
    :return:
    """
    # the order in the list is consistent with the order of appearance in KT file
    p_arms = []
    q_arms = []

    with open(KT_file) as fp_read:
        pass


def Segment_Length(start, end):
    return abs(start - end + 1)
