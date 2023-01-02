# BESO
**BESO** is an efficient implementation for Bi-directional Evolutionary Structural Optimization (BESO) method using C# and C++. The 2D and 3D programs are calibrated with [1] and [4], respectively.


## Dependencies
Noted that all math operations are run on C++.

- [Eigen](https://gitlab.com/libeigen/eigen)
- [Intel® oneAPI Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html#gs.ilgw70)

- If you want to use CHOLMOD (Default), please install [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html) (Windows users can access [suitesparse-metis-for-windows](https://github.com/jlblancoc/suitesparse-metis-for-windows))



## Usage

### 2D BESO
```C#
BESO2D beso = new BESO2D(3.0, 0.5); // filter radius = 3, target volume = 50%
beso.Initialize(80, 50); // resolution 80 x 50

Console.WriteLine(beso.ModelInfo());

// Run optimization
while (!beso.convergence)
{
	beso.Optimize();
	Console.WriteLine(beso.info);
}

Console.ReadKey();
```

Output:

```
=================== Model Info ===================
Disp Nodes: 4131
Disp Elements: 4000

Opt Nodes: 4131
Opt Elements: 4000

=================== Parameters Info ===================
xCount: 80
yCount: 50

Iter: 1, Volume: 0.98, Compliance: 11.573984662017612, Change: 1
Iter: 2, Volume: 0.9603999999999999, Compliance: 11.575623556495295, Change: 1
Iter: 3, Volume: 0.9411919999999999, Compliance: 11.593911784562104, Change: 1
......
```

- If you want to visualize it, please manually output design variables by access `beso.Xe`.

### 3D BESO
Serial mode
```C#
static void Main(string[] args)
{
    Stopwatch sw = new Stopwatch();
    sw.Start();

    BESO3D beso = new BESO3D(3.0, 0.5);
    beso.Initialize(40, 20, 30);
    Console.WriteLine(beso.ModelInfo());

    sw.Stop();
    Console.WriteLine(
        "======================= Init. time: "
        + sw.ElapsedMilliseconds.ToString()
        + " =======================");

    while (!beso.convergence)
    {
        sw.Restart();
        beso.Optimize();
        sw.Stop();
        Console.WriteLine(beso.info);
        Console.WriteLine(
            "======================= It. time: " 
            + sw.ElapsedMilliseconds.ToString() 
            + " =======================");
    }
    Console.ReadKey();
}
```

Output:
```
=================== Model Info ===================
Nodes: 26691
Elements: 24000
Parallel mode: False

=================== Parameters Info ===================
xCount: 40
yCount: 30
zCount: 20

======================= Init. time: 405 =======================
Iter: 1, Volume: 0.97, Compliance: 0.984, Change: 1
======================= It. time: 12480 =======================
Iter: 2, Volume: 0.941, Compliance: 0.984, Change: 1
======================= It. time: 18989 =======================
Iter: 3, Volume: 0.913, Compliance: 0.986, Change: 1
======================= It. time: 19053 =======================
......
```

Parallel mode:
```C#
static void Main(string[] args)
{
    Stopwatch sw = new Stopwatch();
    sw.Start();

    BESO3D beso = new BESO3D(3.0, 0.5);
    beso.Initialize(40, 20, 30, true);
    Console.WriteLine(beso.ModelInfo());

    sw.Stop();
    Console.WriteLine(
        "======================= Init. time: "
        + sw.ElapsedMilliseconds.ToString()
        + " =======================");

    while (!beso.convergence)
    {
        sw.Restart();
        beso.Optimize();
        sw.Stop();
        Console.WriteLine(beso.info);
        Console.WriteLine(
            "======================= It. time: " 
            + sw.ElapsedMilliseconds.ToString() 
            + " =======================");
    }
    Console.ReadKey();
}
```

Output:
```
=================== Model Info ===================
Nodes: 26691
Elements: 24000
Parallel mode: True

=================== Parameters Info ===================
xCount: 40
yCount: 30
zCount: 20

======================= Init. time: 738 =======================
Iter: 1, Volume: 0.97, Compliance: 0.984, Change: 1
======================= It. time: 5516 =======================
Iter: 2, Volume: 0.941, Compliance: 0.984, Change: 1
======================= It. time: 5229 =======================
Iter: 3, Volume: 0.913, Compliance: 0.986, Change: 1
======================= It. time: 5481 =======================
......
```

## References
[1] [CISM_BESO_2D](https://www.cism.org.au/tools)
[2] Zuo, Z.H. and Xie, Y.M., 2015. A simple and compact Python code for complex 3D topology optimization. Advances in Engineering Software, 85, pp.1-11.
[3] Huang, X. and Xie, M., 2010. Evolutionary topology optimization of continuum structures: methods and applications. John Wiley & Sons.
[4] Huang, R. and Huang, X., 2011. Matlab implementation of 3D topology optimization using BESO. Incorporating Sustainable Practice in Mechanics of Structures and Materials, pp.813-818.
[5] Ferrari, F. and Sigmund, O., 2020. A new generation 99 line Matlab code for compliance topology optimization and its extension to 3D. Structural and Multidisciplinary Optimization, 62(4), pp.2211-2228.
