import random


def Segment_Length(start, end):
    return abs(start - end + 1)


def deletion(genome, sequence_start, sequence_end, centromere_start, centromere_end):
    upper_len = centromere_start - sequence_start
    lower_len = sequence_end - centromere_end

    # get length of event
    event_length = random.randint(1000000, 2000000)

    # upper/lower chromosome event: 0 upper, 1 lower
    chr_segment = random.choices([0, 1], weights=[upper_len, lower_len], k=1)

    # find exact sequence location for the event
    if chr_segment == 0:
        event_start_location = random.randrange(sequence_start, centromere_start - event_length)
    else:
        event_start_location = random.randrange(centromere_end + 1, sequence_end - event_length + 1)

    new_genome = genome[:event_start_location] + genome[event_start_location + event_length:]
    return new_genome


