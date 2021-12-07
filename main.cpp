#include <cmath>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <cstdlib> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <cstdio>  //        Remove "-fopenmp" for g++ version < 4.2

#include <array>
#include <vector>
#include <cassert>

// Usage: time ./cleanpt 5000 && display image.ppm

struct Vec {
    double x, y, z; // position, also color (r,g,b)

    Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec operator+(const Vec& b) const {
        return {x + b.x, y + b.y, z + b.z};
    }

    Vec operator-(const Vec& b) const {
        return {x - b.x, y - b.y, z - b.z};
    }

    Vec operator*(double b) const {
        return {x * b, y * b, z * b};
    }

    Vec mult(const Vec& b) const {
        return {x * b.x, y * b.y, z * b.z};
    }

    Vec& norm() {
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }

    double dot(const Vec& b) const {
        return x * b.x + y * b.y + z * b.z;
    }

    // cross:
    Vec operator%(Vec& b) const {
        return {y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x};
    }
};

struct Ray {
    Vec o, d;

    Ray(Vec o_, Vec d_) : o(o_), d(d_) {
    }
};

enum ReflT { DIFF, SPEC, REFR }; // material types, used in radiance()

struct Sphere {
    double rad;  // radius
    Vec p, e, c; // position, emission, color
    ReflT refl;  // reflection type (DIFFuse, SPECular, REFRactive)

    Sphere(double rad_, Vec p_, Vec e_, Vec c_, ReflT refl_)
        : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {
    }

    double intersect(const Ray& r) const { // returns distance, 0 if nohit
        Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t;
        double eps = 1e-4;
        double b = op.dot(r.d);
        double det = b * b - op.dot(op) + rad * rad;
        if (det < 0) {
            return 0;
        } else {
            det = sqrt(det);
        }
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

std::array spheres{
    // Scene: radius, position, emission, color, material
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25),
           DIFF), // Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75),
           DIFF),                                                     // Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF), // Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),       // Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF), // Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75),
           DIFF),                                                      // Top
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC), // Mirr
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR), // Glas
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF) // Lite
};

inline double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int to_int(double x) {
    return static_cast<int>(std::lround(pow(clamp(x), 1 / 2.2) * 255));
}

inline bool intersect(const Ray& r, double& t, std::size_t& id) {
    std::size_t n = spheres.size();
    assert(n > 0);
    double d;
    double inf = t = 1e20;
    for (std::size_t i = n; i != 0; i--) {
        const auto sphere_idx = i - 1;
        if (((d = spheres[sphere_idx].intersect(r)) != 0.0) && d < t) {
            t = d;
            id = sphere_idx;
        }
    }
    return t < inf;
}

Vec radiance(const Ray& r, int depth, unsigned short* Xi) {
    double t;           // distance to intersection
    std::size_t id = 0; // id of intersected object
    if (!intersect(r, t, id)) {
        return {}; // if miss, return black
    }
    const Sphere& obj = spheres[id]; // the hit object
    Vec x = r.o + r.d * t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    Vec f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
    if (++depth > 5) {
        if (erand48(Xi) < p) {
            f = f * (1 / p);
        } else {
            return obj.e; // R.R.
        }
    }
    if (obj.refl == DIFF) { // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi);
        double r2 = erand48(Xi);
        double r2s = sqrt(r2);
        Vec w = nl;
        Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
        Vec v = w % u;
        Vec d =
            (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    } else if (obj.refl == SPEC) // Ideal SPECULAR reflection
        return obj.e +
               f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    Ray refl_ray(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;                 // Ray from outside going in?
    double nc = 1;
    double nt = 1.5;
    double nnt = into ? nc / nt : nt / nc;
    double ddn = r.d.dot(nl);
    double cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) <
        0) // Total internal reflection
        return obj.e + f.mult(radiance(refl_ray, depth, Xi));
    Vec tdir =
        (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc;
    double b = nt + nc;
    double R0 = a * a / (b * b);
    double c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c;
    double Tr = 1 - Re;
    double P = .25 + .5 * Re;
    double RP = Re / P;
    double TP = Tr / (1 - P);
    return obj.e + f.mult(depth > 2
                              ? (erand48(Xi) < P
                                     ? // Russian roulette
                                     radiance(refl_ray, depth, Xi) * RP
                                     : radiance(Ray(x, tdir), depth, Xi) * TP)
                              : radiance(refl_ray, depth, Xi) * Re +
                                    radiance(Ray(x, tdir), depth, Xi) * Tr);
}

int main(int argc, char* argv[]) {
    int w = 1024;
    int h = 768;
    int samps = argc == 2 ? atoi(argv[1]) / 4 : 1;             // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
    Vec cx = Vec(w * .5135 / h);
    Vec cy = (cx % cam.d).norm() * .5135;
    Vec r;
    std::vector<Vec> c;
    c.resize(static_cast<std::size_t>(w) * static_cast<std::size_t>(h));
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
    for (int y = 0; y < h; y++) {                        // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4,
                100. * y / (h - 1));
        unsigned short Xi[3] = {0, 0, static_cast<unsigned short>(y * y * y)};
        for (int x = 0; x < w; x++) { // Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2;
                 sy++) {                                    // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec()) { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++) {
                        double r1 = 2 * erand48(Xi);
                        double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi);
                        double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) +
                                cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) *
                                    (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[static_cast<std::size_t>(i)] =
                        c[static_cast<std::size_t>(i)] +
                        Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
            }
        }
    }
    FILE* f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (std::size_t i = 0;
         i < static_cast<std::size_t>(w) * static_cast<std::size_t>(h); i++) {
        fprintf(f, "%d %d %d ", to_int(c[i].x), to_int(c[i].y), to_int(c[i].z));
    }
}
