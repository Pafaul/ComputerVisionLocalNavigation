from filters import mf
import numpy as np
from PIL import Image
from sys import argv
from uuid import uuid4
import time


def main(path_to_image: str = "test.jpg", sliding_window_size: int = 3) -> None:
    """
    Точка входа в программу.
    :param path_to_image: путь до обрабатываемого изображения.
    :param sliding_window_size: размер скользящего окна.
    :return: сохраняет обработанное изображение.
    """

    img = Image.open(path_to_image).convert("L")
    arr = np.array(img)
    removed_noise = mf(arr, sliding_window_size)
    img = Image.fromarray(removed_noise).convert("L")
    save_name = str(uuid4().hex) + ".png"
    img.save(save_name)


if __name__ == "__main__":
    exception_msg = "Программа запущена с параметрами по умалочанию."
    try:
        path_to_image = argv[1]
        sliding_window_size = int(argv[2])
        print(f"Программа запущенна с параметрами: {path_to_image}, {sliding_window_size}")
        start_time = time.time()
        main(path_to_image, sliding_window_size)
        print(f"time: {(time.time() - start_time)}")
    except FileNotFoundError as e:
        print(exception_msg + " Переданный файл изображения не найден.")
        main()
    except IndexError as e:
        print(exception_msg + " Переданно неверное количество агрументов.")
        main()
