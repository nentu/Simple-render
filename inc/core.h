#ifndef CORE_H
#define CORE_H

#include <embree3/rtcore.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <limits>
#include <chrono>
#include <string>
#include <sstream>

// ==================== Математические структуры ====================

struct Vec3
{
    float x, y, z;

    Vec3();
    Vec3(float x, float y, float z);

    Vec3 operator+(const Vec3 &v) const;
    Vec3 operator-(const Vec3 &v) const;
    Vec3 operator*(float s) const;
    Vec3 operator*(const Vec3 &v) const;
    Vec3 operator/(float s) const;
    Vec3 operator-() const;
    Vec3 &operator+=(const Vec3 &v);
    Vec3 &operator*=(float s);

    float dot(const Vec3 &v) const;
    Vec3 cross(const Vec3 &v) const;
    float length() const;
    float lengthSquared() const;
    Vec3 normalized() const;
    float max_component() const;

    static Vec3 min(const Vec3 &a, const Vec3 &b);
    static Vec3 max(const Vec3 &a, const Vec3 &b);
};

// RGB цвет (используется как для яркости, так и для коэффициентов)
using Color = Vec3;

// ==================== Генератор случайных чисел ====================

class RandomGenerator
{
private:
    std::mt19937 gen;
    std::uniform_real_distribution<float> dist;

    // Преобразование из локальной системы координат (z = нормаль) в мировую
    Vec3 localToWorld(const Vec3 &localDir, const Vec3 &normal);

public:
    RandomGenerator();
    RandomGenerator(unsigned seed);

    float uniform();

    // Равномерная точка на полусфере
    Vec3 uniformHemisphere(const Vec3 &normal);

    // Косинусно-взвешенное направление (для выборки по значимости Ламберта)
    Vec3 cosineHemisphere(const Vec3 &normal);

    // Равномерная точка внутри треугольника
    Vec3 uniformTriangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);
};

// ==================== Материалы ====================

enum class MaterialType
{
    DIFFUSE,  // Диффузное отражение (Ламберт)
    SPECULAR, // Зеркальное отражение
    MIXED,    // Смешанный (диффузный + зеркальный)
    EMISSIVE  // Излучающий (источник света)
};

struct Material
{
    MaterialType type;
    Color diffuseColor;  // Коэффициент диффузного отражения
    Color specularColor; // Коэффициент зеркального отражения
    Color emission;      // Излучение (для источников света)
    float specularProb;  // Вероятность зеркального отражения для смешанного материала

    Material();

    // Создание диффузного материала
    static Material Diffuse(const Color &color);

    // Создание зеркального материала
    static Material Specular(const Color &color = Color(1, 1, 1));

    // Создание смешанного материала
    static Material Mixed(const Color &diffuse, const Color &specular, float specProb);

    // Создание излучающего материала
    static Material Emissive(const Color &emission, const Color &diffuse = Color(0, 0, 0));

    // Проверка физичности: сумма отражений <= 1 для каждого канала
    bool isPhysical() const;
};

// ==================== Треугольник ====================

struct Triangle
{
    Vec3 v0, v1, v2; // Вершины
    Vec3 normal;     // Нормаль
    int materialId;  // Индекс материала
    float area;      // Площадь

    Triangle();
    Triangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, int matId = 0);

    // Получить случайную точку на треугольнике
    Vec3 samplePoint(RandomGenerator &rng) const;
};

// ==================== Источник света ====================

struct AreaLight
{
    Triangle triangle;
    Color emission;
    float power; // Общая мощность источника

    AreaLight(const Triangle &tri, const Color &emit);
};

// ==================== Камера ====================

struct Camera
{
    Vec3 position;
    Vec3 forward;
    Vec3 right;
    Vec3 up;
    float fov; // Угол поля зрения (в радианах)
    int width, height;

    Camera(const Vec3 &pos, const Vec3 &lookAt, const Vec3 &worldUp,
           float fovDegrees, int w, int h);

    // Генерация луча для пикселя (x, y) с антиалиасингом
    void generateRay(int px, int py, float &ox, float &oy, float &oz,
                     float &dx, float &dy, float &dz, RandomGenerator &rng) const;
};

// ==================== Сцена с Embree ====================

class Scene
{
private:
    RTCDevice device;
    RTCScene scene;
    std::vector<Triangle> triangles;
    std::vector<Material> materials;
    std::vector<AreaLight> lights;
    float totalLightPower;

public:
    Scene();
    ~Scene();

    // Добавить материал
    int addMaterial(const Material &mat);

    // Добавить треугольник
    void addTriangle(const Triangle &tri);

    // Добавить четырехугольник (два треугольника)
    void addQuad(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, int matId);

    // Добавить источник света
    void addLight(const Triangle &tri, const Color &emission);

    // Построить BVH после добавления всей геометрии
    void commit();

    // Трассировка луча
    bool intersect(float ox, float oy, float oz,
                   float dx, float dy, float dz,
                   float &t, int &triId) const;

    // Проверка видимости (теневой луч)
    bool isVisible(const Vec3 &from, const Vec3 &to) const;

    const Triangle &getTriangle(int id) const;
    const Material &getMaterial(int id) const;
    const std::vector<AreaLight> &getLights() const;
    float getTotalLightPower() const;
    size_t getTriangleCount() const;

    // Загрузка OBJ файла
    bool loadOBJ(const std::string &filename, int materialId,
                 const Vec3 &offset = Vec3(0, 0, 0), float scale = 1.0f);
};

#endif // CORE_H