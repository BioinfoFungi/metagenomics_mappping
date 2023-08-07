# import tensorflow as tf
# from tensorflow.python.client import device_lib
#
# print("TensorFlow version:", tf.__version__)
#
# print(device_lib.list_local_devices())

import tensorflow as tf
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
print("TensorFlow version:", tf.__version__)

tf.test.is_built_with_cuda()
tf.test.is_gpu_available()
