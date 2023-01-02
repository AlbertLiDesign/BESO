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
Iter: 4, Volume: 0.9223681599999999, Compliance: 11.638274260804778, Change: 1
Iter: 5, Volume: 0.9039207967999998, Compliance: 11.720364729863821, Change: 1
Iter: 6, Volume: 0.8858423808639998, Compliance: 11.829552785502962, Change: 1
Iter: 7, Volume: 0.8681255332467198, Compliance: 12.1207173730191, Change: 1
Iter: 8, Volume: 0.8507630225817854, Compliance: 12.141048701409142, Change: 1
Iter: 9, Volume: 0.8337477621301497, Compliance: 12.266337255276941, Change: 1
Iter: 10, Volume: 0.8170728068875467, Compliance: 12.433127047486394, Change: 1
Iter: 11, Volume: 0.8007313507497957, Compliance: 12.619066548896832, Change: 0.05522096161123532
Iter: 12, Volume: 0.7847167237347998, Compliance: 12.761172749349265, Change: 0.05632890438978203
Iter: 13, Volume: 0.7690223892601038, Compliance: 12.990086968597417, Change: 0.06088873482642055
Iter: 14, Volume: 0.7536419414749017, Compliance: 13.273188062454635, Change: 0.06655712813882725
Iter: 15, Volume: 0.7385691026454037, Compliance: 13.646780870623429, Change: 0.07401635253135297
Iter: 16, Volume: 0.7237977205924956, Compliance: 13.606511431682002, Change: 0.07628159315724707
Iter: 17, Volume: 0.7093217661806457, Compliance: 13.822192980453435, Change: 0.08225564336664914
Iter: 18, Volume: 0.6951353308570327, Compliance: 14.00427890325594, Change: 0.08376691330585334
Iter: 19, Volume: 0.6812326242398921, Compliance: 14.220447965448477, Change: 0.08152067059762082
Iter: 20, Volume: 0.6676079717550942, Compliance: 14.446040430011402, Change: 0.07365836677876518
Iter: 21, Volume: 0.6542558123199923, Compliance: 14.68505948976579, Change: 0.07393552767661239
Iter: 22, Volume: 0.6411706960735924, Compliance: 14.895105093445641, Change: 0.07294716364282203
Iter: 23, Volume: 0.6283472821521205, Compliance: 15.103697017242103, Change: 0.07311165916108432
Iter: 24, Volume: 0.6157803365090782, Compliance: 15.339042353032065, Change: 0.07458465236350045
Iter: 25, Volume: 0.6034647297788965, Compliance: 15.613187719398386, Change: 0.07898233505767902
Iter: 26, Volume: 0.5913954351833186, Compliance: 15.857784843937374, Change: 0.07910865287345392
Iter: 27, Volume: 0.5795675264796523, Compliance: 16.155166311905187, Change: 0.08052417058226405
Iter: 28, Volume: 0.5679761759500592, Compliance: 16.366105036129305, Change: 0.08153930102340529
Iter: 29, Volume: 0.5566166524310581, Compliance: 16.684961905951454, Change: 0.08336712014949008
Iter: 30, Volume: 0.5454843193824369, Compliance: 16.999584701215152, Change: 0.08497941900610889
Iter: 31, Volume: 0.5345746329947881, Compliance: 17.308245306001172, Change: 0.08729787143818145
Iter: 32, Volume: 0.5238831403348924, Compliance: 17.597255030248157, Change: 0.08822047772187364
Iter: 33, Volume: 0.5134054775281945, Compliance: 17.891309887448273, Change: 0.09012926555396449
Iter: 34, Volume: 0.5031373679776306, Compliance: 18.25047163875412, Change: 0.09134749613208464
Iter: 35, Volume: 0.5, Compliance: 18.566561308957834, Change: 0.09200473918688763
Iter: 36, Volume: 0.5, Compliance: 18.655974424614403, Change: 0.08917670555110416
Iter: 37, Volume: 0.5, Compliance: 18.716669592679995, Change: 0.08386484918272692
Iter: 38, Volume: 0.5, Compliance: 18.696091948921552, Change: 0.0740554070582994
Iter: 39, Volume: 0.5, Compliance: 18.676114575494704, Change: 0.05979253427712605
Iter: 40, Volume: 0.5, Compliance: 18.65747860252943, Change: 0.04227567793944596
Iter: 41, Volume: 0.5, Compliance: 18.634828591925494, Change: 0.026600364974053773
Iter: 42, Volume: 0.5, Compliance: 18.657526250493238, Change: 0.013477843356505027
Iter: 43, Volume: 0.5, Compliance: 18.677605486952473, Change: 0.004497832104448301
Iter: 44, Volume: 0.5, Compliance: 18.70753045609539, Change: 0.00025246148204533366
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
yCount: 20
zCount: 30

======================= Init. time: 439 =======================
Iter: 1, Volume: 0.97, Compliance: 1.897, Change: 1
======================= It. time: 12118 =======================
Iter: 2, Volume: 0.941, Compliance: 1.897, Change: 1
======================= It. time: 17979 =======================
Iter: 3, Volume: 0.913, Compliance: 1.897, Change: 1
======================= It. time: 18082 =======================
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
yCount: 20
zCount: 30

======================= Init. time: 725 =======================
Iter: 1, Volume: 0.97, Compliance: 1.897, Change: 1
======================= It. time: 4932 =======================
Iter: 2, Volume: 0.941, Compliance: 1.762, Change: 1
======================= It. time: 4826 =======================
Iter: 3, Volume: 0.913, Compliance: 1.748, Change: 1
======================= It. time: 4783 =======================
......
```

## References

[1] [CISM_BESO_2D](https://www.cism.org.au/tools)
[2] Zuo, Z.H. and Xie, Y.M., 2015. A simple and compact Python code for complex 3D topology optimization. Advances in Engineering Software, 85, pp.1-11.
[3] Huang, X. and Xie, M., 2010. Evolutionary topology optimization of continuum structures: methods and applications. John Wiley & Sons.
[4] Huang, R. and Huang, X., 2011. Matlab implementation of 3D topology optimization using BESO. Incorporating Sustainable Practice in Mechanics of Structures and Materials, pp.813-818.
[5] Ferrari, F. and Sigmund, O., 2020. A new generation 99 line Matlab code for compliance topology optimization and its extension to 3D. Structural and Multidisciplinary Optimization, 62(4), pp.2211-2228.
