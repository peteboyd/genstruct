"""
sample.py

A BTC secondary building unit for sampling purposes
"""

atom_labels = [
          ["C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "C",
           "O",
           "O",
           "O",
           "O",
           "O",
           "O",
           "O",
           "O",
           "H",
           "H",
           "H",
           "H",
           "H",
           "H"],
          ["Zn",
            "X",
            "X",
            "X",
            "X",
            "X"],
          ["C",
           "H",
           "C",
           "H",
           "C",
           "C",
           "H",
           "C",
           "H",
           "C",
           "C",
           "O",
           "O",
           "C",
           "O",
           "O"],
          ["C",  
           "C",  
           "C",  
           "C",  
           "C",  
           "C",  
           "C",  
           "C",  
           "C",  
           "C",  
           "H",  
           "H",  
           "H",  
           "H",  
           "C",  
           "O",  
           "O",  
           "C",  
           "O",  
           "O",  
           "C",  
           "O",  
           "O",  
           "C",  
           "O",  
           "O"],
          ["X",
           "P",
           "O",
           "O",
           "O",
           "O",
           "X"]
            ]

atom_coordinates = [
           [[-2.394,  2.658, -1.359],
           [-1.781,  1.403, -0.701],
           [-0.281,  1.104, -0.873],
           [ 0.602,  2.041, -1.714],
           [-0.009,  3.280, -2.392],
           [-1.508,  3.596, -2.206],
           [-3.893,  2.972, -1.171],
           [-4.770,  2.045, -0.302],
           [-4.159,  0.786,  0.334],
           [-2.673,  0.453,  0.126],
           [-2.119,  4.851, -2.865],
           [-1.226,  5.803, -3.690],
           [ 0.263,  5.470, -3.890],
           [ 0.872,  4.208, -3.256],
           [-3.622,  5.146, -2.697],
           [-4.505,  4.210, -1.855],
           [-6.257,  2.363, -0.045],
           [-2.102, -0.837,  0.751],
           [-1.795,  7.097, -4.314],
           [ 2.361,  3.895, -3.507],
           [-7.089,  1.419,  0.791],
           [-6.883,  3.618, -0.599],
           [-0.661, -1.229,  0.539],
           [-3.004, -1.724,  1.576],
           [ 3.187,  4.848, -4.338],
           [ 2.999,  2.645, -2.952],
           [-0.894,  7.989, -5.135],
           [-3.236,  7.490, -4.099],
           [ 0.182,  0.192, -0.378],
           [ 1.711,  1.812, -1.830],
           [-4.812,  0.095,  0.964],
           [ 0.919,  6.164, -4.514],
           [-4.083,  6.056, -3.194],
           [-5.615,  4.436, -1.741]
           ],
           [[8.213,  6.699,  7.161],
           [ 7.658,  8.007,  8.575],
           [10.147,  6.466,  7.643],
           [ 7.375,  4.925,  7.581],
           [ 8.879,  6.318,  5.308],
           [ 6.735,  7.408,  6.005]
           ],
           [[ 0.830,  1.389,  1.545],
            [ 0.259,  2.306,  1.539],
            [ 2.183,  1.408,  1.903],
            [ 2.656,  2.340,  2.174],
            [ 2.922,  0.219,  1.910],
            [ 2.308, -0.989,  1.559],
            [ 2.878, -1.906,  1.564],
            [ 0.955, -1.008,  1.200],
            [ 0.481, -1.940,  0.929],
            [ 0.216,  0.181,  1.193],
            [-1.273,  0.160,  0.799],
            [-1.869, -1.013,  0.457],
            [-1.991,  1.315,  0.792], 
            [ 4.410,  0.240,  2.304], 
            [ 5.128, -0.915,  2.311], 
            [ 5.007,  1.414,  2.646] 
            ],
           [[-2.690,  1.256, -0.667],
            [-4.221,  1.195, -0.511],
            [-4.858,  0.150,  0.424],
            [-3.964, -0.834,  1.202],
            [-2.433, -0.773,  1.045],
            [-1.796,  0.272,  0.110],
            [-2.052,  2.301, -1.602],
            [-2.946,  3.285, -2.380],
            [-4.477,  3.224, -2.223],
            [-5.114,  2.179, -1.289],
            [-4.436, -1.607,  1.894],
            [-0.663,  0.317, -0.006],
            [-2.474,  4.058, -3.072],
            [-6.247,  2.134, -1.172],
            [-1.540, -1.756,  1.823],
            [-0.188, -1.703,  1.685],
            [-2.103, -2.679,  2.649],
            [-6.389,  0.089,  0.581],
            [-6.951, -0.833,  1.407],
            [-7.178,  0.958, -0.106],
            [-0.522,  2.362, -1.759],
            [ 0.267,  1.493, -1.072],
            [ 0.041,  3.285, -2.585],
            [-5.370,  4.208, -3.001],
            [-6.722,  4.154, -2.863],
            [-4.808,  5.131, -3.827]
            ],
           [[-0.793,  0.628,  1.000],
            [-0.034,  0.640, -1.002],
            [ 1.416,  0.661, -0.618],
            [-0.445, -0.764, -1.331],
            [-0.929,  1.188,  0.233],
            [-0.264,  1.587, -2.297],
            [ 0.265,  1.262, -3.029]]
           ]
