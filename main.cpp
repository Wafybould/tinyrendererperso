#include "tgaimage.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <climits>
#include <limits>
#include "matrix.h"
#include "model.h"
#include "geometry.h"

constexpr int width = 1024;
constexpr int height = 1024;
Model* model = NULL;

void line(int x0, int y0, int x1, int y1, TGAImage &img)
{
    bool steep = std::abs(x0 - x1) < std::abs(y0 - y1);
    if (steep)
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
    }
    if (x1 < x0)
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1-x0;
    int dy = y1-y0;
    int derror = 2*std::abs(dy);
    int error = 0;
    int y = y0;
    for (int x = x0; x < x1; x++)
    {
        if (steep)
        {
            img.set(y, x, TGAColor(255, 255, 255));
        }
        else
        {
            img.set(x, y, TGAColor(255, 255, 255));
        }
        error += derror;
        if (error > dx){
            y += y1 > y0 ? 1 : -1;
            error -= 2*dx;
        }
    }
}

void triangle(Vec3f a, Vec3f b, Vec3f c, TGAImage &framebuffer, TGAColor color, float* zbuffer) {

    //Pas de triangle plat
    if (a.y == b.y && a.y == c.y) return;
    if (a.x == b.x && a.x == c.x) return;

    Vec2f bbmin = {std::numeric_limits<float>::max(),  std::numeric_limits<float>::max()};
    Vec2f bbmax = { -std::numeric_limits<float>::min(), -std::numeric_limits<float>::min() };
    bbmin.x = std::min(bbmin.x, a.x);
    bbmin.x = std::min(bbmin.x, b.x);
    bbmin.x = std::min(bbmin.x, c.x);
    bbmin.y = std::min(bbmin.y, a.y);
    bbmin.y = std::min(bbmin.y, b.y);
    bbmin.y = std::min(bbmin.y, c.y);

    bbmax.x = std::max(bbmax.x, a.x);
    bbmax.x = std::max(bbmax.x, b.x);
    bbmax.x = std::max(bbmax.x, c.x);
    bbmax.y = std::max(bbmax.y, a.y);
    bbmax.y = std::max(bbmax.y, b.y);
    bbmax.y = std::max(bbmax.y, c.y);
    Vec3f P;
    /*matrix33 pts(a.x, b.x, c.x, a.y, b.y, c.y, 1, 1, 1);
    matrix33 invPts = pts.inverse();*/
    for (P.x = bbmin.x; P.x <= bbmax.x; P.x++)
    {
        for (P.y = bbmin.y; P.y <= bbmax.y; P.y++)
        {
            /*std::vector<float> uvw = invPts.mult3x1(std::vector<float>{P.x, P.y, 1 });
            if (uvw[0] < 0 || uvw[1] < 0 || uvw[2] < 0) {
                //std::cout << uvw[2] << std::endl;
                continue;
            }*/
            Vec3f s[2];
            for (int i = 2; i--;) {
                s[i][0] = c[i] - a[i];
                s[i][1] = b[i] - a[i];
                s[i][2] = a[i] - P[i];
            }
            Vec3f u(s[0].y * s[1].z - s[0].z * s[1].y, s[0].z * s[1].x - s[0].x * s[1].z, s[0].x * s[1].y - s[0].y * s[1].x);
            if (std::abs(u.z) > 0) {
                Vec3f res(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
                if (res.x < 0 || res.y < 0 || res.z < 0) continue;
                P.z = a.z * res.x;
                P.z += b.z * res.y;
                P.z += c.z * res.z;
                if (zbuffer[int(P.x + P.y * width)] < P.z) {
                    zbuffer[int(P.x + P.y * width)] = P.z;
                    framebuffer.set(P.x, P.y, color);
                }
            }
        }
    }
}

int main([[maybe_unused]]int argc, [[maybe_unused]]char const *argv[])
{
    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("obj/african_head/african_head.obj");
    }
    /*std::vector<float> pts;
    std::vector<int> tri;
    {
        std::ifstream in;
        in.open(argv[1], std::ifstream::in);
        if (in.fail()) {
            std::cerr << "Failed to open " << argv[1] << std::endl;
            return -1;
        }
        std::string line;
        char trash;
        while (!in.eof()) {
            std::getline(in, line);
            std::istringstream iss(line.c_str());
            if (!line.compare(0, 2, "v ")) {
                iss >> trash;
                for (int i : {0, 1, 2}) {
                    float v;
                    iss >> v;
                    pts.push_back(v);
                }
            }
            if (!line.compare(0, 2, "f ")) {
                iss >> trash;
                int f, t, n;
                while (iss >> f >> trash >> t >> trash >> n) {
                    tri.push_back(f-1);
                }
            }
        }
        in.close();
    }

    int nverts = pts.size() / 3;
    int ntri = tri.size() / 3;*/
    
    /*matrix33 m(1., 2., -1., 2., 1., 2., -1., 2., 1.);
    m.print();
    matrix33 m2 = m.inverse();
    m2.print();*/

    TGAImage framebuffer(width, height, TGAImage::RGB);
    //triangle(Vec2i(0, 1), Vec2i(0, 1024), Vec2i(1023, 1024), framebuffer);
    //triangle(Vec2i(1, 0), Vec2i(1024, 0), Vec2i(1024, 1023), framebuffer);
    
    float* zbuffer = new float[width * height];
    for (int i = width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
    Vec3f light_dir(0, 0, -1);
    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3];
        Vec3f world_coords[3];
        for (int j = 0; j < 3; j++) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec3f(int((v.x + 1.) * width / 2. + .5), int((v.y + 1.) * height / 2. + .5), v.z);
            world_coords[j] = v;
        }
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        float intensity = n * light_dir;
        if (intensity > 0) {
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], framebuffer, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255), zbuffer);
        }
    }

    /*float* zbuffer = new float[width * height];
    for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    TGAImage image(width, height, TGAImage::RGB);
    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec3f pts[3];
        for (int j = 0; j < 3; j++) pts[j] = world2screen(model->vert(face[j]));
        triangle(pts[0], pts[1], pts[2], image, TGAColor(rand() % 255, rand() % 255, rand() % 255, 255), zbuffer);
    }


    /*for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        for (int j = 0; j < 3; j++) {
            Vec3f world_coords = model->vert(face[j]);
            screen_coords[j] = Vec2i((world_coords.x + 1.) * width / 2., (world_coords.y + 1.) * height / 2.);
        }
        triangle(screen_coords[0], screen_coords[1], screen_coords[2], framebuffer);
    }*/

    //line(10, 13, 99, 99, framebuffer);
    /*for (int t = 0; t < ntri; t++)
    {
        for (int s : { 0,1,2 })
        {
            int u = tri[t * 3 + s];
            int v = tri[t * 3 + (s + 1) % 3];
            int x0 = (pts[u * 3 + 0] + 1) / 2 * width;
            int y0 = (pts[u * 3 + 1] + 1) / 2 * height;
            int x1 = (pts[v * 3 + 0] + 1) / 2 * width;
            int y1 = (pts[v * 3 + 1] + 1) / 2 * height;
            line(x0, y0, x1, y1, framebuffer);
        }
    }*/
    framebuffer.write_tga_file("framebuffer.tga");
    return 0;
}
