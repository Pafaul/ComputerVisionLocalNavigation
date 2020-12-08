import numpy as np


def mf(data: np.ndarray, filter_size: int) -> np.ndarray:
    """

    :param data: Преобразованное в numpy массив обрабатываемое изображение.
    :param filter_size: Размер скользящего окна.
    :return: Преобразованное в numpy массив обработанное изображение.
    """

    temp = []
    indexer = filter_size // 2
    data_final = np.zeros((len(data),len(data[0])))
    for i in range(len(data)):

        for j in range(len(data[0])):

            for z in range(filter_size):
                if i + z - indexer < 0 or i + z - indexer > len(data) - 1:
                    for c in range(filter_size):
                        temp.append(0)
                else:
                    if j + z - indexer < 0 or j + indexer > len(data[0]) - 1:
                        temp.append(0)
                    else:
                        for k in range(filter_size):
                            temp.append(data[i + z - indexer][j + k - indexer])

            temp.sort()
            data_final[i][j] = temp[len(temp) // 2]
            temp = []
    return data_final