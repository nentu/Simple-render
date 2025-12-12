#include "scenes.h"
#include <cmath>

void createCylinderAndTetrahedron(Scene &scene)
{
    // Материалы
    int floorMat = scene.addMaterial(Material::Diffuse(Color(0.8f, 0.8f, 0.8f)));
    int cylinderMat = scene.addMaterial(Material::Diffuse(Color(0.7f, 0.3f, 0.2f)));
    int tetraMat = scene.addMaterial(Material::Specular(Color(0.9f, 0.9f, 0.9f)));
    int lightMat = scene.addMaterial(Material::Emissive(Color(12.0f, 12.0f, 10.0f)));

    float size = 4.0f;

    // Пол
    scene.addQuad(
        Vec3(-size, 0, -size), Vec3(size, 0, -size),
        Vec3(size, 0, size), Vec3(-size, 0, size),
        floorMat);

    // Стены
    scene.addQuad(
        Vec3(-size, 0, -size), Vec3(-size, size * 2, -size),
        Vec3(size, size * 2, -size), Vec3(size, 0, -size),
        floorMat);

    // Цилиндр (аппроксимация 8 гранями)
    float cylHeight = 3.0f;
    float cylRadius = 1.0f;
    float cx = 1.5f, cz = 1.0f;
    int segments = 8;

    for (int i = 0; i < segments; i++)
    {
        float angle1 = (i * 2 * M_PI) / segments;
        float angle2 = ((i + 1) * 2 * M_PI) / segments;

        Vec3 v1(cx + cylRadius * cos(angle1), 0, cz + cylRadius * sin(angle1));
        Vec3 v2(cx + cylRadius * cos(angle2), 0, cz + cylRadius * sin(angle2));
        Vec3 v3(cx + cylRadius * cos(angle2), cylHeight, cz + cylRadius * sin(angle2));
        Vec3 v4(cx + cylRadius * cos(angle1), cylHeight, cz + cylRadius * sin(angle1));

        // Боковая грань
        scene.addTriangle(Triangle(v1, v2, v3, cylinderMat));
        scene.addTriangle(Triangle(v1, v3, v4, cylinderMat));

        // Верхняя грань (крышка)
        scene.addTriangle(Triangle(Vec3(cx, cylHeight, cz), v4, v3, cylinderMat));
    }

    // Тетраэдр (пирамида с треугольным основанием)
    float tx = -1.5f, tz = -1.0f, th = 2.5f;
    Vec3 base1(tx - 1.0f, 0, tz - 0.5f);
    Vec3 base2(tx + 1.0f, 0, tz - 0.5f);
    Vec3 base3(tx, 0, tz + 1.0f);
    Vec3 top(tx, th, tz);

    // Грани тетраэдра
    scene.addTriangle(Triangle(base1, base2, top, tetraMat));
    scene.addTriangle(Triangle(base2, base3, top, tetraMat));
    scene.addTriangle(Triangle(base3, base1, top, tetraMat));
    scene.addTriangle(Triangle(base1, base3, base2, tetraMat));

    // Точечный источник света (маленький треугольник)
    Triangle lightTri(
        Vec3(0, size * 2 - 0.1f, 0),
        Vec3(0.2f, size * 2 - 0.1f, 0.2f),
        Vec3(-0.2f, size * 2 - 0.1f, 0.2f),
        lightMat);

    scene.addTriangle(lightTri);
    scene.addLight(lightTri, Color(12.0f, 12.0f, 10.0f));
}

void createChairAndTableScene(Scene &scene)
{
    // Материалы
    int floorMat = scene.addMaterial(Material::Diffuse(Color(0.7f, 0.7f, 0.6f)));
    int tableMat = scene.addMaterial(Material::Diffuse(Color(0.4f, 0.3f, 0.2f)));
    int chairMat = scene.addMaterial(Material::Diffuse(Color(0.3f, 0.2f, 0.1f)));
    int metalMat = scene.addMaterial(Material::Specular(Color(0.8f, 0.8f, 0.8f)));
    int lightMat = scene.addMaterial(Material::Emissive(Color(15.0f, 15.0f, 12.0f)));

    float roomSize = 8.0f;

    // Пол
    scene.addQuad(
        Vec3(-roomSize, 0, -roomSize), Vec3(roomSize, 0, -roomSize),
        Vec3(roomSize, 0, roomSize), Vec3(-roomSize, 0, roomSize),
        floorMat);

    // Стены
    scene.addQuad(
        Vec3(-roomSize, 0, -roomSize), Vec3(-roomSize, roomSize, -roomSize),
        Vec3(roomSize, roomSize, -roomSize), Vec3(roomSize, 0, -roomSize),
        floorMat);

    // Стол
    float tableHeight = 2.5f;
    float tableTopSize = 3.0f;
    float legThickness = 0.2f;

    // Столешница
    scene.addQuad(
        Vec3(-tableTopSize / 2, tableHeight, -tableTopSize / 2),
        Vec3(tableTopSize / 2, tableHeight, -tableTopSize / 2),
        Vec3(tableTopSize / 2, tableHeight, tableTopSize / 2),
        Vec3(-tableTopSize / 2, tableHeight, tableTopSize / 2),
        tableMat);

    // Ножки стола
    scene.addQuad(
        Vec3(-tableTopSize / 2 + legThickness, 0, -tableTopSize / 2 + legThickness),
        Vec3(-tableTopSize / 2 + legThickness, tableHeight, -tableTopSize / 2 + legThickness),
        Vec3(-tableTopSize / 2, tableHeight, -tableTopSize / 2 + legThickness),
        Vec3(-tableTopSize / 2, 0, -tableTopSize / 2 + legThickness),
        tableMat);

    scene.addQuad(
        Vec3(tableTopSize / 2 - legThickness, 0, -tableTopSize / 2 + legThickness),
        Vec3(tableTopSize / 2 - legThickness, tableHeight, -tableTopSize / 2 + legThickness),
        Vec3(tableTopSize / 2, tableHeight, -tableTopSize / 2 + legThickness),
        Vec3(tableTopSize / 2, 0, -tableTopSize / 2 + legThickness),
        tableMat);

    scene.addQuad(
        Vec3(-tableTopSize / 2 + legThickness, 0, tableTopSize / 2 - legThickness),
        Vec3(-tableTopSize / 2 + legThickness, tableHeight, tableTopSize / 2 - legThickness),
        Vec3(-tableTopSize / 2, tableHeight, tableTopSize / 2 - legThickness),
        Vec3(-tableTopSize / 2, 0, tableTopSize / 2 - legThickness),
        tableMat);

    scene.addQuad(
        Vec3(tableTopSize / 2 - legThickness, 0, tableTopSize / 2 - legThickness),
        Vec3(tableTopSize / 2 - legThickness, tableHeight, tableTopSize / 2 - legThickness),
        Vec3(tableTopSize / 2, tableHeight, tableTopSize / 2 - legThickness),
        Vec3(tableTopSize / 2, 0, tableTopSize / 2 - legThickness),
        tableMat);

    // Стул
    float chairHeight = 1.5f;
    float chairSeatHeight = 1.0f;
    float chairSeatSize = 1.2f;
    float chairBackHeight = 2.0f;

    // Сиденье стула
    scene.addQuad(
        Vec3(-chairSeatSize / 2, chairSeatHeight, -chairSeatSize / 2),
        Vec3(chairSeatSize / 2, chairSeatHeight, -chairSeatSize / 2),
        Vec3(chairSeatSize / 2, chairSeatHeight, chairSeatSize / 2),
        Vec3(-chairSeatSize / 2, chairSeatHeight, chairSeatSize / 2),
        chairMat);

    // Спинка стула
    scene.addQuad(
        Vec3(-chairSeatSize / 2, chairSeatHeight, chairSeatSize / 2),
        Vec3(chairSeatSize / 2, chairSeatHeight, chairSeatSize / 2),
        Vec3(chairSeatSize / 2, chairBackHeight, chairSeatSize / 2),
        Vec3(-chairSeatSize / 2, chairBackHeight, chairSeatSize / 2),
        chairMat);

    // Ножки стула (металлические)
    scene.addQuad(
        Vec3(-chairSeatSize / 2 + 0.1f, 0, -chairSeatSize / 2 + 0.1f),
        Vec3(-chairSeatSize / 2 + 0.1f, chairSeatHeight, -chairSeatSize / 2 + 0.1f),
        Vec3(-chairSeatSize / 2 + 0.2f, chairSeatHeight, -chairSeatSize / 2 + 0.1f),
        Vec3(-chairSeatSize / 2 + 0.2f, 0, -chairSeatSize / 2 + 0.1f),
        metalMat);

    scene.addQuad(
        Vec3(chairSeatSize / 2 - 0.2f, 0, -chairSeatSize / 2 + 0.1f),
        Vec3(chairSeatSize / 2 - 0.2f, chairSeatHeight, -chairSeatSize / 2 + 0.1f),
        Vec3(chairSeatSize / 2 - 0.1f, chairSeatHeight, -chairSeatSize / 2 + 0.1f),
        Vec3(chairSeatSize / 2 - 0.1f, 0, -chairSeatSize / 2 + 0.1f),
        metalMat);

    scene.addQuad(
        Vec3(-chairSeatSize / 2 + 0.1f, 0, chairSeatSize / 2 - 0.2f),
        Vec3(-chairSeatSize / 2 + 0.1f, chairSeatHeight, chairSeatSize / 2 - 0.2f),
        Vec3(-chairSeatSize / 2 + 0.2f, chairSeatHeight, chairSeatSize / 2 - 0.2f),
        Vec3(-chairSeatSize / 2 + 0.2f, 0, chairSeatSize / 2 - 0.2f),
        metalMat);

    // Источник света на потолке
    float lightSize = 1.5f;
    float lightY = roomSize - 0.1f;

    Triangle lightTri1(
        Vec3(-lightSize, lightY, -lightSize),
        Vec3(lightSize, lightY, -lightSize),
        Vec3(lightSize, lightY, lightSize),
        lightMat);
    Triangle lightTri2(
        Vec3(-lightSize, lightY, -lightSize),
        Vec3(lightSize, lightY, lightSize),
        Vec3(-lightSize, lightY, lightSize),
        lightMat);

    scene.addTriangle(lightTri1);
    scene.addTriangle(lightTri2);
    scene.addLight(lightTri1, Color(15.0f, 15.0f, 12.0f));
    scene.addLight(lightTri2, Color(15.0f, 15.0f, 12.0f));
}