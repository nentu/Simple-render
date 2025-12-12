/**
 * Лабораторная работа №5: Трассировка путей (Path Tracing)
 *
 * Цель: Освоить простейшие методы синтеза изображений трехмерных сцен
 * с учетом глобального освещения методом трассировки путей.
 *
 * Реализовано:
 * - Трассировка путей с глобальным освещением
 * - Геометрия сцены на основе треугольных сеток (Embree API)
 * - Материалы: диффузное (Ламберт) и зеркальное отражение
 * - Выборка по значимости и русская рулетка
 * - Протяженные источники света с диаграммой Ламберта
 * - Точечная камера с антиалиасингом
 * - HDR -> LDR преобразование с гамма-коррекцией
 * - Вывод в формате PPM
 */

#include "pathtracer.h"
#include "utils.h"
#include "scenes.h"
#include <iostream>
#include <string>

// ==================== Вспомогательные функции ====================

void printUsage(const char *programName)
{
    std::cout << "Использование: " << programName << " [опции]" << std::endl;
    std::cout << "Опции:" << std::endl;
    std::cout << "  -w <ширина>      Ширина изображения (по умолчанию: 512)" << std::endl;
    std::cout << "  -h <высота>      Высота изображения (по умолчанию: 512)" << std::endl;
    std::cout << "  -s <samples>     Число лучей на пиксель (по умолчанию: 64)" << std::endl;
    std::cout << "  -d <depth>       Максимальная глубина трассировки (по умолчанию: 10)" << std::endl;
    std::cout << "  -o <файл>        Имя выходного файла (по умолчанию: output.ppm)" << std::endl;
    std::cout << "  -e <exposure>    Экспозиция (по умолчанию: авто)" << std::endl;
    std::cout << "  -g <gamma>       Гамма-коррекция (по умолчанию: 2.2)" << std::endl;
    std::cout << "  --obj <файл>     Загрузить OBJ модель" << std::endl;
    std::cout << "  --scale <число>  Масштаб OBJ модели (по умолчанию: 1.0)" << std::endl;
    std::cout << "  --scene <номер>  Номер сцены (1-6, по умолчанию: 4)" << std::endl;
    std::cout << "  --help           Показать эту справку" << std::endl;
    std::cout << std::endl;
    std::cout << "Доступные сцены:" << std::endl;
    std::cout << "  1 - Cornell Box (классическая)" << std::endl;
    std::cout << "  2 - Mirror Room (комната с зеркалами)" << std::endl;
    std::cout << "  3 - Three Spheres (три сферы)" << std::endl;
    std::cout << "  4 - Colored Lights (цветное освещение)" << std::endl;
    std::cout << "  5 - Cylinder and Tetrahedron (цилиндр и тетраэдр)" << std::endl;
    std::cout << "  6 - Chair and Table (стул и стол)" << std::endl;
}

// ==================== Main ====================

int main(int argc, char *argv[])
{
    // Параметры по умолчанию
    int width = 512;
    int height = 512;
    int samplesPerPixel = 64;
    int maxDepth = 10;
    int sceneNumber = 4; // Цветное освещение по умолчанию
    std::string outputFile = "output.ppm";
    std::string objFile = "";
    float objScale = 1.0f;
    float exposure = -1.0f; // -1 означает автоэкспозицию
    float gamma = 2.2f;

    // Парсинг аргументов командной строки
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "--help")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if (arg == "-w" && i + 1 < argc)
        {
            width = std::atoi(argv[++i]);
        }
        else if (arg == "-h" && i + 1 < argc)
        {
            height = std::atoi(argv[++i]);
        }
        else if (arg == "-s" && i + 1 < argc)
        {
            samplesPerPixel = std::atoi(argv[++i]);
        }
        else if (arg == "-d" && i + 1 < argc)
        {
            maxDepth = std::atoi(argv[++i]);
        }
        else if (arg == "-o" && i + 1 < argc)
        {
            outputFile = argv[++i];
        }
        else if (arg == "-e" && i + 1 < argc)
        {
            exposure = std::atof(argv[++i]);
        }
        else if (arg == "-g" && i + 1 < argc)
        {
            gamma = std::atof(argv[++i]);
        }
        else if (arg == "--obj" && i + 1 < argc)
        {
            objFile = argv[++i];
        }
        else if (arg == "--scale" && i + 1 < argc)
        {
            objScale = std::atof(argv[++i]);
        }
        else if (arg == "--scene" && i + 1 < argc)
        {
            sceneNumber = std::atoi(argv[++i]);
        }
    }

    // Проверка параметров
    width = std::max(100, std::min(width, 2000));
    height = std::max(100, std::min(height, 2000));
    samplesPerPixel = std::max(1, std::min(samplesPerPixel, 10000));
    maxDepth = std::max(1, std::min(maxDepth, 100));
    sceneNumber = std::max(1, std::min(sceneNumber, 6));

    std::cout << "========================================" << std::endl;
    std::cout << "     Path Tracer — Работа №5" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Настройки рендеринга:" << std::endl;
    std::cout << "  Размер изображения: " << width << " x " << height << std::endl;
    std::cout << "  Сэмплов на пиксель: " << samplesPerPixel << std::endl;
    std::cout << "  Глубина трассировки: " << maxDepth << std::endl;
    std::cout << "  Файл результата: " << outputFile << std::endl;
    std::cout << "  Коррекция гаммы: " << gamma << std::endl;
    std::cout << "  Выбранная сцена: " << sceneNumber << std::endl;
    if (!objFile.empty())
    {
        std::cout << "  OBJ файл: " << objFile << " (масштаб: " << objScale << ")" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    // Создание сцены
    Scene scene;

    if (!objFile.empty())
    {
        // Создаём простую сцену с OBJ моделью
        int whiteMat = scene.addMaterial(Material::Diffuse(Color(0.73f, 0.73f, 0.73f)));
        int objMat = scene.addMaterial(Material::Diffuse(Color(0.8f, 0.3f, 0.3f)));
        int lightMat = scene.addMaterial(Material::Emissive(Color(15.0f, 15.0f, 15.0f)));

        // Пол
        float floorSize = 10.0f;
        scene.addQuad(
            Vec3(-floorSize, 0, -floorSize), Vec3(floorSize, 0, -floorSize),
            Vec3(floorSize, 0, floorSize), Vec3(-floorSize, 0, floorSize),
            whiteMat);

        // Загружаем OBJ модель
        scene.loadOBJ(objFile, objMat, Vec3(0, 0, 0), objScale);

        // Источник света
        Triangle light1(Vec3(-3, 10, -3), Vec3(3, 10, -3), Vec3(3, 10, 3), lightMat);
        Triangle light2(Vec3(-3, 10, -3), Vec3(3, 10, 3), Vec3(-3, 10, 3), lightMat);
        scene.addTriangle(light1);
        scene.addTriangle(light2);
        scene.addLight(light1, Color(15.0f, 15.0f, 15.0f));
        scene.addLight(light2, Color(15.0f, 15.0f, 15.0f));
    }
    else
    {
        // Выбор сцены по номеру
        switch (sceneNumber)
        {
        case 1:
            createCylinderAndTetrahedron(scene);
            break;
        case 2:
            createChairAndTableScene(scene);
            break;
        default:
            createCylinderAndTetrahedron(scene);
            break;
        }
    }

    scene.commit();

    std::cout << "Сцена создана: " << scene.getTriangleCount() << " треугольников" << std::endl;

    // Создание камеры (разная позиция для OBJ и предопределенных сцен)
    Vec3 cameraPos, cameraLookAt;
    Vec3 cameraUp(0, 1, 0);

    if (!objFile.empty())
    {
        // Для OBJ модели: камера смотрит на центр сцены
        cameraPos = Vec3(5.0f * objScale, 3.0f * objScale, 5.0f * objScale);
        cameraLookAt = Vec3(0, 1.0f * objScale, 0);
    }
    else
    {
        // Для предопределенных сцен
        switch (sceneNumber)
        {
        case 1: // Cornell Box
            cameraPos = Vec3(0, 5.0f, 14.0f);
            cameraLookAt = Vec3(0, 5.0f, 0);
            break;
        case 2: // Mirror Room
            cameraPos = Vec3(0, 5.0f, 12.0f);
            cameraLookAt = Vec3(0, 5.0f, 0);
            break;
        case 3: // Three Spheres
            cameraPos = Vec3(0, 3.0f, 10.0f);
            cameraLookAt = Vec3(0, 2.0f, 0);
            break;
        case 4: // Colored Lights
            cameraPos = Vec3(0, 5.0f, 12.0f);
            cameraLookAt = Vec3(0, 3.0f, 0);
            break;
        case 5: // Cylinder and Tetrahedron
            cameraPos = Vec3(0, 4.0f, 10.0f);
            cameraLookAt = Vec3(0, 2.0f, 0);
            break;
        case 6: // Chair and Table
            cameraPos = Vec3(0, 3.0f, 8.0f);
            cameraLookAt = Vec3(0, 1.5f, 0);
            break;
        default:
            cameraPos = Vec3(0, 5.0f, 12.0f);
            cameraLookAt = Vec3(0, 3.0f, 0);
            break;
        }
    }

    Camera camera(cameraPos, cameraLookAt, cameraUp, 45.0f, width, height);

    // Рендеринг
    PathTracer tracer(scene, camera, maxDepth, samplesPerPixel);
    std::vector<Color> hdrImage = tracer.render();

    // Постобработка и сохранение
    bool autoExposure = (exposure < 0);
    if (exposure < 0)
        exposure = 1.0f;

    std::vector<unsigned char> ldrImage = tonemapAndGammaCorrect(hdrImage, exposure, gamma, autoExposure);
    savePPM(outputFile, ldrImage, width, height);

    std::cout << "========================================" << std::endl;
    std::cout << "    Рендеринг завершен!" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}