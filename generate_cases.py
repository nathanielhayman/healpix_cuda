import numpy as np
import cv2

for i in range(5):
    img = np.asarray([
        np.asarray(255 * np.ones(10**i))
        for x in range(0, 10**i)
    ])

    cv2.imwrite(f"../src/image{10**i}x{10**i}.jpg", img)