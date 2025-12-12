#include "core.h"

// ==================== Реализации Vec3 ====================

Vec3::Vec3() : x(0), y(0), z(0) {}
Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

Vec3 Vec3::operator+(const Vec3 &v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
Vec3 Vec3::operator-(const Vec3 &v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
Vec3 Vec3::operator*(float s) const { return Vec3(x * s, y * s, z * s); }
Vec3 Vec3::operator*(const Vec3 &v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
Vec3 Vec3::operator/(float s) const { return Vec3(x / s, y / s, z / s); }
Vec3 Vec3::operator-() const { return Vec3(-x, -y, -z); }
Vec3 &Vec3::operator+=(const Vec3 &v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}
Vec3 &Vec3::operator*=(float s)
{
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

float Vec3::dot(const Vec3 &v) const { return x * v.x + y * v.y + z * v.z; }

Vec3 Vec3::cross(const Vec3 &v) const
{
    return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}

float Vec3::length() const { return std::sqrt(x * x + y * y + z * z); }
float Vec3::lengthSquared() const { return x * x + y * y + z * z; }

Vec3 Vec3::normalized() const
{
    float len = length();
    if (len > 0)
        return *this / len;
    return *this;
}

float Vec3::max_component() const { return std::max({x, y, z}); }

Vec3 Vec3::min(const Vec3 &a, const Vec3 &b)
{
    return Vec3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

Vec3 Vec3::max(const Vec3 &a, const Vec3 &b)
{
    return Vec3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

// ==================== Реализации RandomGenerator ====================

RandomGenerator::RandomGenerator() : gen(std::random_device{}()), dist(0.0f, 1.0f) {}
RandomGenerator::RandomGenerator(unsigned seed) : gen(seed), dist(0.0f, 1.0f) {}

float RandomGenerator::uniform() { return dist(gen); }

// Равномерная точка на полусфере
Vec3 RandomGenerator::uniformHemisphere(const Vec3 &normal)
{
    float z = uniform();
    float r = std::sqrt(std::max(0.0f, 1.0f - z * z));
    float phi = 2.0f * M_PI * uniform();

    Vec3 localDir(r * std::cos(phi), r * std::sin(phi), z);
    return localToWorld(localDir, normal);
}

// Косинусно-взвешенное направление (для выборки по значимости Ламберта)
Vec3 RandomGenerator::cosineHemisphere(const Vec3 &normal)
{
    float r1 = uniform();
    float r2 = uniform();

    float sinTheta = std::sqrt(r1);
    float cosTheta = std::sqrt(1.0f - r1);
    float phi = 2.0f * M_PI * r2;

    Vec3 localDir(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    return localToWorld(localDir, normal);
}

// Равномерная точка внутри треугольника
Vec3 RandomGenerator::uniformTriangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
{
    float r1 = uniform();
    float r2 = uniform();

    if (r1 + r2 > 1.0f)
    {
        r1 = 1.0f - r1;
        r2 = 1.0f - r2;
    }

    return v0 * (1.0f - r1 - r2) + v1 * r1 + v2 * r2;
}

// Преобразование из локальной системы координат (z = нормаль) в мировую
Vec3 RandomGenerator::localToWorld(const Vec3 &localDir, const Vec3 &normal)
{
    Vec3 n = normal.normalized();
    Vec3 tangent = (std::abs(n.x) > 0.9f) ? Vec3(0, 1, 0) : Vec3(1, 0, 0);
    Vec3 bitangent = n.cross(tangent).normalized();
    tangent = bitangent.cross(n);

    return tangent * localDir.x + bitangent * localDir.y + n * localDir.z;
}

// ==================== Реализации Material ====================

Material::Material() : type(MaterialType::DIFFUSE),
                       diffuseColor(0.8f, 0.8f, 0.8f),
                       specularColor(1.0f, 1.0f, 1.0f),
                       emission(0, 0, 0),
                       specularProb(0.0f) {}

// Создание диффузного материала
Material Material::Diffuse(const Color &color)
{
    Material m;
    m.type = MaterialType::DIFFUSE;
    m.diffuseColor = color;
    return m;
}

// Создание зеркального материала
Material Material::Specular(const Color &color)
{
    Material m;
    m.type = MaterialType::SPECULAR;
    m.specularColor = color;
    return m;
}

// Создание смешанного материала
Material Material::Mixed(const Color &diffuse, const Color &specular, float specProb)
{
    Material m;
    m.type = MaterialType::MIXED;
    m.diffuseColor = diffuse;
    m.specularColor = specular;
    m.specularProb = std::min(std::max(specProb, 0.0f), 1.0f);
    return m;
}

// Создание излучающего материала
Material Material::Emissive(const Color &emission, const Color &diffuse)
{
    Material m;
    m.type = MaterialType::EMISSIVE;
    m.emission = emission;
    m.diffuseColor = diffuse;
    return m;
}

// Проверка физичности: сумма отражений <= 1 для каждого канала
bool Material::isPhysical() const
{
    Color total = diffuseColor + specularColor;
    return total.x <= 1.0f && total.y <= 1.0f && total.z <= 1.0f;
}

// ==================== Реализации Triangle ====================

Triangle::Triangle() : materialId(0), area(0) {}

Triangle::Triangle(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, int matId)
    : v0(v0), v1(v1), v2(v2), materialId(matId)
{
    Vec3 edge1 = v1 - v0;
    Vec3 edge2 = v2 - v0;
    Vec3 cross = edge1.cross(edge2);
    normal = cross.normalized();
    area = cross.length() * 0.5f;
}

// Получить случайную точку на треугольнике
Vec3 Triangle::samplePoint(RandomGenerator &rng) const
{
    return rng.uniformTriangle(v0, v1, v2);
}

// ==================== Реализации AreaLight ====================

AreaLight::AreaLight(const Triangle &tri, const Color &emit)
    : triangle(tri), emission(emit)
{
    // Мощность = интеграл по площади (emission * cos(theta) * dA)
    // Для равномерного Ламберта: power = emission * area * pi
    power = (emit.x + emit.y + emit.z) / 3.0f * tri.area * M_PI;
}

// ==================== Реализации Camera ====================

Camera::Camera(const Vec3 &pos, const Vec3 &lookAt, const Vec3 &worldUp,
               float fovDegrees, int w, int h)
    : position(pos), width(w), height(h)
{
    fov = fovDegrees * M_PI / 180.0f;
    forward = (lookAt - pos).normalized();
    right = forward.cross(worldUp).normalized();
    up = right.cross(forward);
}

// Генерация луча для пикселя (x, y) с антиалиасингом
void Camera::generateRay(int px, int py, float &ox, float &oy, float &oz,
                         float &dx, float &dy, float &dz, RandomGenerator &rng) const
{
    // Случайное смещение внутри пикселя для антиалиасинга
    float jitterX = rng.uniform() - 0.5f;
    float jitterY = rng.uniform() - 0.5f;

    float x = (2.0f * (px + 0.5f + jitterX) / width - 1.0f);
    float y = (1.0f - 2.0f * (py + 0.5f + jitterY) / height);

    float aspectRatio = float(width) / float(height);
    float tanHalfFov = std::tan(fov * 0.5f);

    x *= aspectRatio * tanHalfFov;
    y *= tanHalfFov;

    Vec3 dir = (forward + right * x + up * y).normalized();

    ox = position.x;
    oy = position.y;
    oz = position.z;
    dx = dir.x;
    dy = dir.y;
    dz = dir.z;
}

// ==================== Реализации Scene ====================

// Thread-local контексты для каждого потока
RTCIntersectContext &Scene::getIntersectContext()
{
    static thread_local RTCIntersectContext context;
    static thread_local bool initialized = false;
    if (!initialized)
    {
        rtcInitIntersectContext(&context);
        context.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
        initialized = true;
    }
    return context;
}

RTCIntersectContext &Scene::getShadowContext()
{
    static thread_local RTCIntersectContext context;
    static thread_local bool initialized = false;
    if (!initialized)
    {
        rtcInitIntersectContext(&context);
        context.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
        initialized = true;
    }
    return context;
}

Scene::Scene() : totalLightPower(0)
{
    device = rtcNewDevice(nullptr);
    if (!device)
    {
        std::cerr << "Ошибка: не удалось создать Embree device" << std::endl;
        exit(1);
    }
    scene = rtcNewScene(device);
}

Scene::~Scene()
{
    rtcReleaseScene(scene);
    rtcReleaseDevice(device);
}

// Добавить материал
int Scene::addMaterial(const Material &mat)
{
    materials.push_back(mat);
    return materials.size() - 1;
}

// Добавить треугольник
void Scene::addTriangle(const Triangle &tri)
{
    triangles.push_back(tri);
}

// Добавить четырехугольник (два треугольника)
void Scene::addQuad(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec3 &v3, int matId)
{
    addTriangle(Triangle(v0, v1, v2, matId));
    addTriangle(Triangle(v0, v2, v3, matId));
}

// Добавить источник света
void Scene::addLight(const Triangle &tri, const Color &emission)
{
    AreaLight light(tri, emission);
    lights.push_back(light);
    totalLightPower += light.power;
}

// Построить BVH после добавления всей геометрии
void Scene::commit()
{
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    // Выделяем буферы вершин и индексов
    float *vertices = (float *)rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        3 * sizeof(float), triangles.size() * 3);

    unsigned *indices = (unsigned *)rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
        3 * sizeof(unsigned), triangles.size());

    // Заполняем буферы
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        vertices[i * 9 + 0] = triangles[i].v0.x;
        vertices[i * 9 + 1] = triangles[i].v0.y;
        vertices[i * 9 + 2] = triangles[i].v0.z;
        vertices[i * 9 + 3] = triangles[i].v1.x;
        vertices[i * 9 + 4] = triangles[i].v1.y;
        vertices[i * 9 + 5] = triangles[i].v1.z;
        vertices[i * 9 + 6] = triangles[i].v2.x;
        vertices[i * 9 + 7] = triangles[i].v2.y;
        vertices[i * 9 + 8] = triangles[i].v2.z;

        indices[i * 3 + 0] = i * 3 + 0;
        indices[i * 3 + 1] = i * 3 + 1;
        indices[i * 3 + 2] = i * 3 + 2;
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);
}

// Трассировка луча
bool Scene::intersect(float ox, float oy, float oz,
                      float dx, float dy, float dz,
                      float &t, int &triId) const
{
    RTCRayHit rayhit;
    rayhit.ray.org_x = ox;
    rayhit.ray.org_y = oy;
    rayhit.ray.org_z = oz;
    rayhit.ray.dir_x = dx;
    rayhit.ray.dir_y = dy;
    rayhit.ray.dir_z = dz;
    rayhit.ray.tnear = 0.001f;
    rayhit.ray.tfar = std::numeric_limits<float>::infinity();
    rayhit.ray.mask = -1;
    rayhit.ray.flags = 0;
    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.instID[1] = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.instID[2] = RTC_INVALID_GEOMETRY_ID;

    // Использовать thread_local контекст вместо статического метода
    static thread_local RTCIntersectContext context;
    static thread_local bool contextInitialized = false;
    if (!contextInitialized)
    {
        rtcInitIntersectContext(&context);
        context.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
        contextInitialized = true;
    }

    rtcIntersect1(scene, &context, &rayhit);

    if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
    {
        t = rayhit.ray.tfar;
        triId = rayhit.hit.primID;
        return true;
    }
    return false;
}

// Проверка видимости (теневой луч)
bool Scene::isVisible(const Vec3 &from, const Vec3 &to) const
{
    Vec3 dir = to - from;
    float dist = dir.length();
    dir = dir / dist;

    RTCRay ray;
    ray.org_x = from.x;
    ray.org_y = from.y;
    ray.org_z = from.z;
    ray.dir_x = dir.x;
    ray.dir_y = dir.y;
    ray.dir_z = dir.z;
    ray.tnear = 0.001f;
    ray.tfar = dist - 0.001f;
    ray.mask = -1;
    ray.flags = 0;

    // Использовать тот же контекст что и в intersect
    static thread_local RTCIntersectContext context;
    static thread_local bool contextInitialized = false;
    if (!contextInitialized)
    {
        rtcInitIntersectContext(&context);
        context.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
        contextInitialized = true;
    }

    rtcOccluded1(scene, &context, &ray);

    return ray.tfar >= 0; // Если tfar < 0, луч был заблокирован
}

const Triangle &Scene::getTriangle(int id) const { return triangles[id]; }
const Material &Scene::getMaterial(int id) const { return materials[id]; }
const std::vector<AreaLight> &Scene::getLights() const { return lights; }
float Scene::getTotalLightPower() const { return totalLightPower; }
size_t Scene::getTriangleCount() const { return triangles.size(); }

// Загрузка OBJ файла
bool Scene::loadOBJ(const std::string &filename, int materialId,
                    const Vec3 &offset, float scale)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть OBJ файл: " << filename << std::endl;
        return false;
    }

    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    int triangleCount = 0;

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v")
        {
            // Вершина
            float x, y, z;
            iss >> x >> y >> z;
            vertices.push_back(Vec3(x * scale + offset.x,
                                    y * scale + offset.y,
                                    z * scale + offset.z));
        }
        else if (prefix == "vn")
        {
            // Нормаль
            float nx, ny, nz;
            iss >> nx >> ny >> nz;
            normals.push_back(Vec3(nx, ny, nz));
        }
        else if (prefix == "f")
        {
            // Грань (face) - поддержка форматов: v, v/vt, v/vt/vn, v//vn
            std::vector<int> faceVertices;
            std::string vertexData;

            while (iss >> vertexData)
            {
                int vIdx = 0;
                size_t pos = vertexData.find('/');
                if (pos != std::string::npos)
                {
                    vIdx = std::stoi(vertexData.substr(0, pos));
                }
                else
                {
                    vIdx = std::stoi(vertexData);
                }
                // OBJ индексы начинаются с 1
                faceVertices.push_back(vIdx - 1);
            }

            // Триангуляция полигона (fan triangulation)
            for (size_t i = 1; i + 1 < faceVertices.size(); ++i)
            {
                int i0 = faceVertices[0];
                int i1 = faceVertices[i];
                int i2 = faceVertices[i + 1];

                if (i0 >= 0 && i0 < (int)vertices.size() &&
                    i1 >= 0 && i1 < (int)vertices.size() &&
                    i2 >= 0 && i2 < (int)vertices.size())
                {
                    addTriangle(Triangle(vertices[i0], vertices[i1], vertices[i2], materialId));
                    triangleCount++;
                }
            }
        }
    }

    file.close();
    std::cout << "Загружен OBJ: " << filename << " (" << vertices.size()
              << " вершин, " << triangleCount << " треугольников)" << std::endl;
    return true;
}