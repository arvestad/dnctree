
def words(file):
    for line in file:
        for word in line.strip().split():
            yield word


def read_taxa(file):
    taxa = []
    file_words = words(file)
    taxa_count = int(next(file_words))
    for i in range(taxa_count):
        taxa.append(next(file_words))
    return taxa


def read_distances(file):
    distances = []
    file_words = words(file)
    next(file_words)  # Ignore total rows
    next(file_words)  # Ignore total columns
    populated_rows = int(next(file_words))
    populated_columns = int(next(file_words))
    for i in range(populated_rows * populated_columns):
        distances.append(float(next(file_words)))
    return distances
