**(Last update: 23 Aug, 2014)**

![http://aobench.googlecode.com/hg/ao.png](http://aobench.googlecode.com/hg/ao.png)

[JavaScript+HTML5 version](http://kioku.sys-k.net/aobench_jsx/aobench.html)

[WebGL version(Need Chrome)](http://aobench.googlecode.com/hg/webgl/ao.html)

**aoench** is a small ambient occlusion renderer for benchmarking realworld floating point performance in various languages.

aobench is originally coded in "Proce55ing":http://processing.org/.  aoench consists just about 400 lines of code and requires only standard math functions, therefore it is easy to understand and port aobench for your favorite language.

Follow me at twitter for updates: http://twitter.com/aobench

## Media and publications ##

<a href='http://www.youtube.com/watch?feature=player_embedded&v=FAlfKFdoj7w' target='_blank'><img src='http://img.youtube.com/vi/FAlfKFdoj7w/0.jpg' width='425' height=344 /></a>

[Revision 2012](https://code.google.com/p/aobench/source/detail?r=2012) (PC-64k intro) - SCOTTIE by qatnonoil & kioku/System K & Machia ( aobench was used as a loading screen for demoscene )

<a href='http://www.youtube.com/watch?feature=player_embedded&v=0q8BCBy_9iU' target='_blank'><img src='http://img.youtube.com/vi/0q8BCBy_9iU/0.jpg' width='425' height=344 /></a>

aobench debut at SC11!

  * aobench at TokyoDemoFest 2012
    * http://tokyo-demo-fest.jpn.org/2012/upload/aobench.pdf

<a href='http://www.youtube.com/watch?feature=player_embedded&v=iS1CcTMNHC4' target='_blank'><img src='http://img.youtube.com/vi/iS1CcTMNHC4/0.jpg' width='425' height=344 /></a>

> Extreme aobench at TokyoDemoFest 2013!

  * Bringing SIMD-128 to JavaScript
    * http://esdiscuss.org/notes/2014-07/simd-128-tc39.pdf

Also, aobench is successfully cited in...

  * Extending a C-like Language for Portable SIMD Programming(PPoPP 2012)
    * http://www.cdl.uni-saarland.de/projects/vecimp/vecimp_tr.pdf
  * ispc: A SPMD Compiler for High-Performance CPU Programming, INPAR 2012
    * http://pl887.pairlitesite.com/papers/ispc/ispc_inpar_2012.pdf
  * Whole-Function Vectorization
    * http://www.intel-vci.uni-saarland.de/uploads/tx_sibibtex/10_01.pdf
  * Ubiquitous Virtual Machine Supporting Semantic Types
    * http://www.soumu.go.jp/main_sosiki/joho_tsusin/scope/event/h21yokousyu/session6/wakate6.pdf
  * Improving Memory Management for ruby
    * http://www.atdot.net/~ko1/activities/rubymem.pdf
  * A proposal of stack-based garbage collection and its evaluation in scripting language Lua
    * http://almond.cs.uec.ac.jp/papers/pdf/2010/Komuro_Pro.pdf
  * C++ Language Constructs for Parallel Programming
    * http://www.open-std.org/JTC1/SC22/wg14/www/docs/n1612.pdf
  * Sierra: A SIMD Extension for C++
    * http://www.cdl.uni-saarland.de/papers/lhh14.pdf
  * Your paper here ;-)

## Languages and performances ##


Here's list of languages aoench was ported to. Your effort to porting aobench to another language or reporting performance of aobench is appreciate(let me know your work by following "aobench@twitter":http://twitter.com/aobench or by mail).

  * "Proce55ing":http://lucille.atso-net.jp/blog/?p=632 (Original version) 9.7 sec on Core2 2.16GHz
  * "C":http://lucille.svn.sourceforge.net/viewvc/lucille/angelina/proce55ing/c_reference/ 2.6 sec on Core2 2.16GHz(gcc-4.4 -O3)
    * "VC++ version":AOBenchOnVC9\_090214-0126.zip by wk. 1.45 sec on Core2 3.16GHz.
    * "VC++ version[version](multithreaded.md)":http://d.hatena.ne.jp/Dycoon/20090221 by "Dycoon":http://d.hatena.ne.jp/Dycoon/ . 1.3 sec on Athlon X2 4800+(4 threaded).
  * "NativeClient":http://lucille.atso-net.jp/nacl/ao.html ("source":http://lucille.svn.sourceforge.net/viewvc/lucille/angelina/nacl/ao/ )
  * "Flash10[ActionScript3](ActionScript3.md)":http://d.hatena.ne.jp/keim_at_Si/20090218#p1 by "keim\_at\_Si":http://d.hatena.ne.jp/keim_at_Si . 7.04 sec on Core2 2.16GHz
  * "JavaScript":http://lucille.atso-net.jp/blog/?p=642
  * "Alchemy":http://lucille.atso-net.jp/blog/?p=661 480 sec on Core2 2.16GHz
  * "iPhone":http://lucille.atso-net.jp/blog/?p=656 38 sec
  * "MAXScript":http://lucille.atso-net.jp/blog/?p=718 by "Guillermo":http://www.evvisual.com/ 651 sec
  * "Android DalvikVM":http://d.hatena.ne.jp/takuma104/20081216 by "takuma":http://d.hatena.ne.jp/takuma104/ 2966 sec
  * "Android C":http://katsu-net.blogspot.com/2009/02/android-dev-phone-1-aobench-ni.html by "katsu":http://katsu-net.blogspot.com 229 sec(Note: float precision).
  * "GPUAO(GLSL) in his 4K gfx tool":http://kioku.sys-k.net/archives/2009/01/gpuaoglsl.html by  "kioku":http://kioku.sys-k.net/. 0.01 sec on 8800GT. Ultimetely fast!
  * "MEL":http://blog.taikomatsu.com/2009/02/08/ao-bench-meets-mel/ by "tai":http://blog.taikomatsu.com/. 110 sec on Athlon64 X2 4800+
  * "Ruby":http://d.hatena.ne.jp/miura1729/20090208/1234084920 by "Miura":http://d.hatena.ne.jp/miura1729. 556 sec ("yarv2llvm":http://github.com/miura1729/yarv2llvm/tree/master executes in 200 sec, which is about x3 faster than ruby1.9)
  * "Ypsilon[Scheme](Scheme.md)":http://d.hatena.ne.jp/mjt/20090206/p1 by "mjt":http://d.hatena.ne.jp/mjt
  * "Mosh[Scheme](Scheme.md)":http://d.hatena.ne.jp/higepon/20090208/1234094156 by "higepon":http://d.hatena.ne.jp/higepon
  * "Gauche[Scheme](Scheme.md)":http://mono.kmc.gr.jp/~yhara/d/?date=20090219#p02 by "yhara":http://mono.kmc.gr.jp/~yhara/d . 307 sec on Core2 2.4GHz.
  * "Gauche[Scheme](Scheme.md)":http://practical-scheme.net/wiliki/wiliki.cgi?Gauche:AOBench  by "Shiro":http://practical-scheme.net/wiliki/wiliki.cgi?Shiro . 85 sec on Core2 2.4GHz(26 sec for 4 threaded version).
  * "Delphi":aobench.dpr by "bee":http://www.bee-www.com/. 2.9 sec on Core2 2.4 GHz(almost identital with C version).
  * "C#":AOBenchOnCSConsole\_src\_090213-1611.zip by wk. 3.73 sec on PenM 1.86 GHz.
  * "XNA":http://d.hatena.ne.jp/XELF/20090212 by "XELF":http://d.hatena.ne.jp/XELF
  * "Python":http://mglab.blogspot.com/2009/02/ambient-occlusion-with-python.html by   "rezoo":http://mglab.blogspot.com . 463 sec on Athlon64 X2 3800+.
  * "Python[Fastest](Fastest.md)":http://leonardo-m.livejournal.com/79346.html by "leonardo":http://leonardo-m.livejournal.com . 138.64 sec on Core2 2GHz.
  * "Python+Psyco":http://leonardo-m.livejournal.com/79346.html by leonardo. 29.75 sec on Core2 2GHz(16.72 for 2 threaded version)
  * "Haskell":http://d.hatena.ne.jp/mokehehe/20090214 by "mokehehe":http://d.hatena.ne.jp/mokehehe/ . 21.6 sec on Core2 2.16GHz.
  * "CommonLisp[sbcl](sbcl.md)":http://d.hatena.ne.jp/youz/20090220/1235101933 by "youz":http://d.hatena.ne.jp/youz . 14.62 sec on Core2 1.8GHz.
  * "TJS2":http://d.hatena.ne.jp/shirakusa/20090215/1234712727 by "shirakusa":http://d.hatena.ne.jp/shirakusa . 1324 sec on Pen4 3.0GHz.
  * "lua":http://blog.goo.ne.jp/c5h12/e/8268be99b5a1f821f7e7cbbd7bc1e268 by "c5h12":http://blog.goo.ne.jp/c5h12 . 198 sec on Core2 2.2GHz.
  * "lua[mdiapp](mdiapp.md)":http://d.hatena.ne.jp/ryocotan/20090209/p1 by "ryokotan":http://d.hatena.ne.jp/ryocotan .
  * "Smalltalk":http://d.hatena.ne.jp/sumim/20090227/p1 by "sumim":http://d.hatena.ne.jp/sumim/ . 65 sec on Core2 2.4GHz.
  * "Xtal":http://d.hatena.ne.jp/xtalco/20090227#1235737637 by "xtalco":http://d.hatena.ne.jp/xtalco . 498 sec.
  * "Java[iAppli](iAppli.md)":http://www.pannoki.com/blog/article.php?id=336 by "eitaro":http://www.pannoki.com/blog . 655.36 sec on N-04.
  * "Silverlight2":http://www.geocities.jp/cubic_831/AOBenchOnSilverlight2/index.html by wk. 3.38 sec on Core2 3.16GHz + Firefox3.
  * "F#":http://d.hatena.ne.jp/ototoi/20090308 by "ototoi":http://d.hatena.ne.jp/ototoi/ . 9.4 sec on Core2 2.1GHz.
  * "HSP":http://d.hatena.ne.jp/MATSUZAKI/20090323 by "MATSUZAKI":http://d.hatena.ne.jp/MATSUZAKI . 259 sec on Core2 3GHz.
  * "D":http://leonardo-m.livejournal.com/79346.html by "leonardo":http://leonardo-m.livejournal.com/ . 3.67 sec on Core2 2GHz.
  * "D, using LDC" by "leonardo":http://leonardo-m.livejournal.com/ . 2.95 on Core2 2GHz(faster than C).
  * "Java":http://leonardo-m.livejournal.com/79346.html by leonardo. 6.81 sec on Core2 2GHz.
  * "Java[Parallel Utiltiy](Using.md)":http://hopson.web.fc2.com/aobench/java/ by hopson. 1.68 sec on Core2 2.4GHz(2 threads, float).
  * "OpenCL[SnowLeopard](SnowLeopard.md)":http://kioku.sys-k.net/archives/2009/08/opencl_ao_bench.html by kioku. 2.8 fps on 9400M.
  * "VC++[VS2010's ParallelFor](Using.md)":ao\_vs2010\_parallelfor.cpp by "hiroyuk":http://blogs.msdn.com/hiroyuk/ . 0.642 sec on Core2 Quad Q9000 2.00GHz (VC++ 2010 Beta1, Windows 7)
  * "WebGL":http://aobench.googlecode.com/hg/webgl/ao.html by ando
  * "ISPC":http://ispc.github.com/ by mpp
  * "Ocaml":http://d.hatena.ne.jp/h013/20110625/1309114394 by h013. 6.6 sec on Core2 2.66GHz.
  * "Scala":https://gist.github.com/1029695 http://d.hatena.ne.jp/MATSUZAKI/20110616/1308243608 by MATSUZAKI.  9.4 secs on Core2 Duo P8700 2.53GHz
  * "clang-interpreter":http://d.hatena.ne.jp/ohtorii/20110724 by ohtorii.
  * "Cilk plus": http://software.intel.com/en-us/articles/data-and-thread-parallelism/ by Intel.
  * "js-of-ocaml:http://peppermint.jp/temp/ao/ by zakki.
  * "JSX":http://kioku.sys-k.net/aobench_jsx/ by system-k.
  * "haXe":https://github.com/yoshihiro503/aobench_haxe by yoshihiro503.
  * "PS Vita":http://gyabo.sakura.ne.jp/prog_ao.html by gyabo. 12 secs.
  * "Amiga 1200":http://twitter.com/bpoint42/status/211416863727550464/photo/1 by bpoint42. 36m44s
    * https://github.com/bpoint/aobench-amiga
  * "WebCL":http://www.eaflux.com/aobench_webcl/ by Koji.
  * "Kuin":https://gist.github.com/3482315 http://twitter.com/santarh/status/239973847225495552/photo/1 by santarh
  * "Parallella":https://github.com/parallella/parallella-examples/tree/master/aobench

### C version code ###

```
/*

aobench C code is licensed under 2-clause BSD.

Copyright 2009-2014, Syoyo Fujita
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#define WIDTH        256
#define HEIGHT       256
#define NSUBSAMPLES  2
#define NAO_SAMPLES  8

typedef struct _vec
{
    double x;
    double y;
    double z;
} vec;


typedef struct _Isect
{
    double t;
    vec    p;
    vec    n;
    int    hit; 
} Isect;

typedef struct _Sphere
{
    vec    center;
    double radius;

} Sphere;

typedef struct _Plane
{
    vec    p;
    vec    n;

} Plane;

typedef struct _Ray
{
    vec    org;
    vec    dir;
} Ray;

Sphere spheres[3];
Plane  plane;

static double vdot(vec v0, vec v1)
{
    return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

static void vcross(vec *c, vec v0, vec v1)
{
    
    c->x = v0.y * v1.z - v0.z * v1.y;
    c->y = v0.z * v1.x - v0.x * v1.z;
    c->z = v0.x * v1.y - v0.y * v1.x;
}

static void vnormalize(vec *c)
{
    double length = sqrt(vdot((*c), (*c)));

    if (fabs(length) > 1.0e-17) {
        c->x /= length;
        c->y /= length;
        c->z /= length;
    }
}

void
ray_sphere_intersect(Isect *isect, const Ray *ray, const Sphere *sphere)
{
    vec rs;

    rs.x = ray->org.x - sphere->center.x;
    rs.y = ray->org.y - sphere->center.y;
    rs.z = ray->org.z - sphere->center.z;

    double B = vdot(rs, ray->dir);
    double C = vdot(rs, rs) - sphere->radius * sphere->radius;
    double D = B * B - C;

    if (D > 0.0) {
        double t = -B - sqrt(D);
        
        if ((t > 0.0) && (t < isect->t)) {
            isect->t = t;
            isect->hit = 1;
            
            isect->p.x = ray->org.x + ray->dir.x * t;
            isect->p.y = ray->org.y + ray->dir.y * t;
            isect->p.z = ray->org.z + ray->dir.z * t;

            isect->n.x = isect->p.x - sphere->center.x;
            isect->n.y = isect->p.y - sphere->center.y;
            isect->n.z = isect->p.z - sphere->center.z;

            vnormalize(&(isect->n));
        }
    }
}

void
ray_plane_intersect(Isect *isect, const Ray *ray, const Plane *plane)
{
    double d = -vdot(plane->p, plane->n);
    double v = vdot(ray->dir, plane->n);

    if (fabs(v) < 1.0e-17) return;

    double t = -(vdot(ray->org, plane->n) + d) / v;

    if ((t > 0.0) && (t < isect->t)) {
        isect->t = t;
        isect->hit = 1;
        
        isect->p.x = ray->org.x + ray->dir.x * t;
        isect->p.y = ray->org.y + ray->dir.y * t;
        isect->p.z = ray->org.z + ray->dir.z * t;

        isect->n = plane->n;
    }
}

void
orthoBasis(vec *basis, vec n)
{
    basis[2] = n;
    basis[1].x = 0.0; basis[1].y = 0.0; basis[1].z = 0.0;

    if ((n.x < 0.6) && (n.x > -0.6)) {
        basis[1].x = 1.0;
    } else if ((n.y < 0.6) && (n.y > -0.6)) {
        basis[1].y = 1.0;
    } else if ((n.z < 0.6) && (n.z > -0.6)) {
        basis[1].z = 1.0;
    } else {
        basis[1].x = 1.0;
    }

    vcross(&basis[0], basis[1], basis[2]);
    vnormalize(&basis[0]);

    vcross(&basis[1], basis[2], basis[0]);
    vnormalize(&basis[1]);
}


void ambient_occlusion(vec *col, const Isect *isect)
{
    int    i, j;
    int    ntheta = NAO_SAMPLES;
    int    nphi   = NAO_SAMPLES;
    double eps = 0.0001;

    vec p;

    p.x = isect->p.x + eps * isect->n.x;
    p.y = isect->p.y + eps * isect->n.y;
    p.z = isect->p.z + eps * isect->n.z;

    vec basis[3];
    orthoBasis(basis, isect->n);

    double occlusion = 0.0;

    for (j = 0; j < ntheta; j++) {
        for (i = 0; i < nphi; i++) {
            double theta = sqrt(drand48());
            double phi   = 2.0 * M_PI * drand48();

            double x = cos(phi) * theta;
            double y = sin(phi) * theta;
            double z = sqrt(1.0 - theta * theta);

            // local -> global
            double rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
            double ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
            double rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;

            Ray ray;

            ray.org = p;
            ray.dir.x = rx;
            ray.dir.y = ry;
            ray.dir.z = rz;

            Isect occIsect;
            occIsect.t   = 1.0e+17;
            occIsect.hit = 0;

            ray_sphere_intersect(&occIsect, &ray, &spheres[0]); 
            ray_sphere_intersect(&occIsect, &ray, &spheres[1]); 
            ray_sphere_intersect(&occIsect, &ray, &spheres[2]); 
            ray_plane_intersect (&occIsect, &ray, &plane); 

            if (occIsect.hit) occlusion += 1.0;
            
        }
    }

    occlusion = (ntheta * nphi - occlusion) / (double)(ntheta * nphi);

    col->x = occlusion;
    col->y = occlusion;
    col->z = occlusion;
}

unsigned char
clamp(double f)
{
  int i = (int)(f * 255.5);

  if (i < 0) i = 0;
  if (i > 255) i = 255;

  return (unsigned char)i;
}


void
render(unsigned char *img, int w, int h, int nsubsamples)
{
    int x, y;
    int u, v;

    double *fimg = (double *)malloc(sizeof(double) * w * h * 3);
    memset((void *)fimg, 0, sizeof(double) * w * h * 3);

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            
            for (v = 0; v < nsubsamples; v++) {
                for (u = 0; u < nsubsamples; u++) {
                    double px = (x + (u / (double)nsubsamples) - (w / 2.0)) / (w / 2.0);
                    double py = -(y + (v / (double)nsubsamples) - (h / 2.0)) / (h / 2.0);

                    Ray ray;

                    ray.org.x = 0.0;
                    ray.org.y = 0.0;
                    ray.org.z = 0.0;

                    ray.dir.x = px;
                    ray.dir.y = py;
                    ray.dir.z = -1.0;
                    vnormalize(&(ray.dir));

                    Isect isect;
                    isect.t   = 1.0e+17;
                    isect.hit = 0;

                    ray_sphere_intersect(&isect, &ray, &spheres[0]);
                    ray_sphere_intersect(&isect, &ray, &spheres[1]);
                    ray_sphere_intersect(&isect, &ray, &spheres[2]);
                    ray_plane_intersect (&isect, &ray, &plane);

                    if (isect.hit) {
                        vec col;
                        ambient_occlusion(&col, &isect);

                        fimg[3 * (y * w + x) + 0] += col.x;
                        fimg[3 * (y * w + x) + 1] += col.y;
                        fimg[3 * (y * w + x) + 2] += col.z;
                    }

                }
            }

            fimg[3 * (y * w + x) + 0] /= (double)(nsubsamples * nsubsamples);
            fimg[3 * (y * w + x) + 1] /= (double)(nsubsamples * nsubsamples);
            fimg[3 * (y * w + x) + 2] /= (double)(nsubsamples * nsubsamples);
        
            img[3 * (y * w + x) + 0] = clamp(fimg[3 *(y * w + x) + 0]);
            img[3 * (y * w + x) + 1] = clamp(fimg[3 *(y * w + x) + 1]);
            img[3 * (y * w + x) + 2] = clamp(fimg[3 *(y * w + x) + 2]);
        }
    }

}

void
init_scene()
{
    spheres[0].center.x = -2.0;
    spheres[0].center.y =  0.0;
    spheres[0].center.z = -3.5;
    spheres[0].radius = 0.5;
    
    spheres[1].center.x = -0.5;
    spheres[1].center.y =  0.0;
    spheres[1].center.z = -3.0;
    spheres[1].radius = 0.5;
    
    spheres[2].center.x =  1.0;
    spheres[2].center.y =  0.0;
    spheres[2].center.z = -2.2;
    spheres[2].radius = 0.5;

    plane.p.x = 0.0;
    plane.p.y = -0.5;
    plane.p.z = 0.0;

    plane.n.x = 0.0;
    plane.n.y = 1.0;
    plane.n.z = 0.0;

}

void
saveppm(const char *fname, int w, int h, unsigned char *img)
{
    FILE *fp;

    fp = fopen(fname, "wb");
    assert(fp);

    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", w, h);
    fprintf(fp, "255\n");
    fwrite(img, w * h * 3, 1, fp);
    fclose(fp);
}

int
main(int argc, char **argv)
{
    unsigned char *img = (unsigned char *)malloc(WIDTH * HEIGHT * 3);

    init_scene();

    render(img, WIDTH, HEIGHT, NSUBSAMPLES);

    saveppm("ao.ppm", WIDTH, HEIGHT, img); 

    return 0;
}
```

### Business card aobench ###

by @h013

```
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
typedef double D;typedef int I;struct V
{V operator +(V r){return V(x+r.x,y+r.y
,z+r.z);}V operator*(D r){return V(x*r,
y*r,z*r);}D operator%(V r){return x*r.x
+y*r.y+z*r.z;}D x,y,z;V(){}V operator^(
V r){return V(y*r.z-z*r.y,z*r.x-x*r.z,x
*r.y-y*r.x);}V operator!(){return*this*
(1/sqrt(*this%*this));}V(D a,D b,D c){x
=a;y=b;z=c;}};I W=256,H=256,Q=2,A=8;V S
[3]={V(-2,0,-3.5),V(-0.5,0,-3),V(1,0,-
2.2)};D R[3]={0.5,0.5,0.5};V P(0,-0.5,0
),N(0,1,0);I T(V o,V d,V&p,V&n){D t=1e9
;D b=d%N,c=(o%N+P%N)/b;if(fabs(b)>1e-17
&&0<c){t=c;p=o+d*t;n=N;}for(I i=0;i<3;
++i){V s=o+(S[i]*-1);c=s%d;c=-c-sqrt(c*
c-s%s+R[i]*R[i]);if(0<c&&c<t){t=c;p=o+d
*t;n=!(p+S[i]*-1);}}return t<1e9?1:0;}D
F(){return(D)rand()/RAND_MAX;}I main(){
printf("P3 %d %d 255 ",W,H);D x,y,u,v,X
,Y,e,t,s,i,j,c;t=s=A;for(y=0;y<H;++y)
for(x=0;c=0,x<W;printf("%d %d %d ",(I)c
,(I)c,(I)c),x++)for(v=0;v<Q;++v)for(u=0
;u<Q;++u){V o(0,0,0),d(2*(x+u/Q)/W-1,1-
2*(y+v/Q)/H,-1),p,n;if(T(o,!d,p,n)){e=t
*s;for(j=0;j<t;++j)for(i=0;i<s;X=sqrt(F
()),Y=6.2831*F(),o=!(V(1,0,0)^n),d=!V((
o*cos(Y)*X+(n^o)*sin(Y)*X+n*sqrt(1-X*X)
)),T(p,d,o,o)&&e--,++i);c+=e*255/(Q*Q*t
*s);}}}
```

## Images ##

![http://kioku.sys-k.net/images/gpuao-thumb.jpg](http://kioku.sys-k.net/images/gpuao-thumb.jpg)
http://kioku.sys-k.net/images/gpuao.html

![http://f.hatena.ne.jp/images/fotolife/X/XELF/20090212/20090212115649.png](http://f.hatena.ne.jp/images/fotolife/X/XELF/20090212/20090212115649.png)
http://f.hatena.ne.jp/XELF/20090212115649

http://www.pannoki.com/blog/resources/AObench-Iappli.JPG

![http://hopson.web.fc2.com/aobench/java/images/aobench_mac_java_f_p.png](http://hopson.web.fc2.com/aobench/java/images/aobench_mac_java_f_p.png)

![http://kioku.sys-k.net/images/aocl_cpu.png](http://kioku.sys-k.net/images/aocl_cpu.png)

![http://kioku.sys-k.net/images/aocl_gpu.png](http://kioku.sys-k.net/images/aocl_gpu.png)

![https://p.twimg.com/Au8ahHaCAAAJ4kf.jpg](https://p.twimg.com/Au8ahHaCAAAJ4kf.jpg)