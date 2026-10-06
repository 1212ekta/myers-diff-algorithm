import sys
import gc
from array import array


def read_lines(path):
    """Read a file as raw bytes and split it into lines (keeps any \\r)."""
    with open(path, "rb") as f:
        data = f.read()
    lines = data.split(b"\n")
    # A final newline leaves an empty last piece, so drop it.
    # An empty file becomes [] this way too.
    if lines and lines[-1] == b"":
        lines.pop()
    return lines


def myers_core(seq_a, seq_b, index_shift=0):
    """
    Myers O(ND) diff on two sequences (lists of ints, or strings).
    Returns the edit script as a list of (operation, index):
      ('=', i)  keep seq_a[i]
      ('-', i)  delete seq_a[i]
      ('+', j)  insert seq_b[j]
    index_shift is added to every index in the result. diff() uses it to turn
    positions inside the middle part back into positions in the full sequences.
    """
    len_a, len_b = len(seq_a), len(seq_b)
    max_edits = len_a + len_b
    offset = max_edits + 1                  # +1 leaves room for the sentinel slot on the left
    # furthest_x[diagonal + offset] = furthest x reached on that diagonal
    # (diagonal = x - y). All zeros at the start.
    furthest_x = [0] * (2 * max_edits + 3)
    # saved_rounds[e] = the furthest x of diagonals -e, -e+2, ..., +e right
    # after round e (e+1 numbers, stored compactly). The backward pass needs
    # these to retrace the path.
    saved_rounds = []
    total_edits = 0
    reached_end = False

    # ---- forward pass: try 0 edits, then 1 edit, then 2 edits, ... ----
    for edit_count in range(max_edits + 1):
        # Sentinels: put -1 just outside the diagonals of this round. Then the
        # first diagonal (-edit_count) always chooses "insert" and the last one
        # (+edit_count) always chooses "delete", with no special case in the loop.
        furthest_x[offset - edit_count - 1] = -1
        furthest_x[offset + edit_count + 1] = -1

        # With exactly edit_count edits we can only be on diagonals
        # -edit_count, -edit_count+2, ..., +edit_count.
        for diagonal in range(-edit_count, edit_count + 1, 2):
            slot = offset + diagonal
            # Choose where we came from: the neighbour diagonal that got further.
            if furthest_x[slot - 1] < furthest_x[slot + 1]:
                x = furthest_x[slot + 1]            # step down  = insert from B
            else:
                x = furthest_x[slot - 1] + 1        # step right = delete from A
            y = x - diagonal

            # Follow the "snake": slide along the diagonal while lines match (free).
            while x < len_a and y < len_b and seq_a[x] == seq_b[y]:
                x += 1
                y += 1

            furthest_x[slot] = x

            if x >= len_a and y >= len_b:       # reached the bottom-right corner
                total_edits = edit_count
                reached_end = True
                break
        if reached_end:
            break
        # Save this round for the backward pass: every second slot, i.e. exactly
        # the diagonals -e, -e+2, ..., +e (the last round is not saved).
        saved_rounds.append(
            array("i", furthest_x[offset - edit_count: offset + edit_count + 1: 2])
        )

    # ---- backward pass: walk from the end to the start using saved_rounds ----
    reversed_script = []
    x, y = len_a, len_b
    for edit_count in range(total_edits, 0, -1):
        previous_round = saved_rounds[edit_count - 1]
        # previous_round holds diagonals -(e-1), -(e-1)+2, ..., (e-1), so
        # diagonal d is stored at position (d + (e-1)) // 2.
        round_shift = edit_count - 1
        diagonal = x - y

        # Same choice as the forward pass: which diagonal did we come from?
        if diagonal == -edit_count or (
            diagonal != edit_count
            and previous_round[(diagonal - 1 + round_shift) // 2]
                < previous_round[(diagonal + 1 + round_shift) // 2]
        ):
            prev_diagonal = diagonal + 1
            came_down = True                # last move was an insert
        else:
            prev_diagonal = diagonal - 1
            came_down = False               # last move was a delete

        prev_x = previous_round[(prev_diagonal + round_shift) // 2]
        prev_y = prev_x - prev_diagonal

        # Point right after the single insert/delete move (where the snake begins).
        if came_down:
            snake_start_x, snake_start_y = prev_x, prev_y + 1
        else:
            snake_start_x, snake_start_y = prev_x + 1, prev_y

        # Walk back along the snake: these lines are kept.
        while x > snake_start_x and y > snake_start_y:
            x -= 1
            y -= 1
            reversed_script.append(('=', x + index_shift))

        # Record the one insert or delete that led here.
        if came_down:
            reversed_script.append(('+', prev_y + index_shift))
        else:
            reversed_script.append(('-', prev_x + index_shift))
        x, y = prev_x, prev_y

    # The snake at the very start (before any edit): all kept lines.
    while x > 0:
        x -= 1
        reversed_script.append(('=', x + index_shift))

    reversed_script.reverse()
    return reversed_script


def diff(seq_a, seq_b):
    """Same as myers_core, but first skips the common start and end (a big speed-up)."""
    len_a, len_b = len(seq_a), len(seq_b)

    # Length of the common beginning.
    prefix_len = 0
    while prefix_len < len_a and prefix_len < len_b and seq_a[prefix_len] == seq_b[prefix_len]:
        prefix_len += 1

    # Length of the common ending (not overlapping the prefix).
    suffix_len = 0
    while (suffix_len < len_a - prefix_len and suffix_len < len_b - prefix_len
           and seq_a[len_a - 1 - suffix_len] == seq_b[len_b - 1 - suffix_len]):
        suffix_len += 1

    edit_script = [('=', i) for i in range(prefix_len)]

    # Run Myers only on the middle part; index_shift puts the indexes back into
    # full-sequence positions.
    middle_a = seq_a[prefix_len:len_a - suffix_len]
    middle_b = seq_b[prefix_len:len_b - suffix_len]
    edit_script.extend(myers_core(middle_a, middle_b, prefix_len))

    edit_script.extend(('=', i) for i in range(len_a - suffix_len, len_a))
    return edit_script


def format_ranges(positions):
    """[3,4,5,9] -> '3-6,9-10' (end is exclusive). Empty list -> '.'"""
    if not positions:
        return "."
    parts = []
    range_start = previous = positions[0]
    for position in positions[1:]:
        if position == previous + 1:        # touching: extend the current range
            previous = position
        else:                               # gap: close the range, start a new one
            parts.append(f"{range_start}-{previous + 1}")
            range_start = previous = position
    parts.append(f"{range_start}-{previous + 1}")
    return ",".join(parts)


def highlight_line(old_line, new_line):
    """Build the '? old_ranges | new_ranges' line for one changed pair."""
    old_text = old_line.decode("utf-8")     # str indexes by code point, as required
    new_text = new_line.decode("utf-8")
    deleted_positions = []                  # character positions removed from old_text
    inserted_positions = []                 # character positions added in new_text
    for operation, index in diff(old_text, new_text):
        if operation == '-':
            deleted_positions.append(index)
        elif operation == '+':
            inserted_positions.append(index)
    line = f"? {format_ranges(deleted_positions)} | {format_ranges(inserted_positions)}"
    return line.encode("utf-8")


def main() -> int:
    gc.disable()        # many tuples/lists are created; skip the repeated cleanup scans
    if len(sys.argv) != 4 or sys.argv[1] not in ("lines", "highlight"):
        print("usage: main.py lines|highlight A_PATH B_PATH", file=sys.stderr)
        return 2
    command, a_path, b_path = sys.argv[1:]

    try:
        old_lines = read_lines(a_path)
        new_lines = read_lines(b_path)
    except OSError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    # Give every distinct line a number so comparing lines is a fast int compare.
    line_ids = {}
    old_ids = [line_ids.setdefault(line, len(line_ids)) for line in old_lines]
    new_ids = [line_ids.setdefault(line, len(line_ids)) for line in new_lines]
    del line_ids        # not needed any more; frees the big dictionary

    edit_script = diff(old_ids, new_ids)
    del old_ids, new_ids

    output = []
    pending_deletes = []    # indexes (into old_lines) of '-' lines in the current change block
    pending_inserts = []    # indexes (into new_lines) of '+' lines in the current change block

    def flush_change_block():
        """Print the block: all '-' first, then all '+' (with '?' line for each pair)."""
        for old_index in pending_deletes:
            output.append(b"-" + old_lines[old_index])
        for pair_number, new_index in enumerate(pending_inserts):
            output.append(b"+" + new_lines[new_index])
            # The n-th '+' pairs with the n-th '-'; leftover '+' lines get no '?' line.
            if command == "highlight" and pair_number < len(pending_deletes):
                old_index = pending_deletes[pair_number]
                output.append(highlight_line(old_lines[old_index], new_lines[new_index]))
        pending_deletes.clear()
        pending_inserts.clear()

    for operation, index in edit_script:
        if operation == '=':
            if pending_deletes or pending_inserts:
                flush_change_block()            # a kept line ends the current change block
            output.append(b" " + old_lines[index])
        elif operation == '-':
            pending_deletes.append(index)
        else:
            pending_inserts.append(index)
    flush_change_block()                    # the file may end inside a change block

    if output:
        sys.stdout.buffer.write(b"\n".join(output) + b"\n")
    return 0


raise SystemExit(main())
