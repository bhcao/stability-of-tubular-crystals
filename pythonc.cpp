/*********************************************************************//**
 * AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)
 * 
 * @file pythonc.cpp
 * @brief Python 库生成文件
 * @version 0.2.0
 * 
 * Python 库生成文件，会使用此文件描述 Python/C API 库的函数
 **************************************************************************/

#include <Python.h>

#include "energy.h"
#include "model.h"
#include "molecule.h"
#include "narray.h"
#include "nmath.h"

//! 暂时只支持单例模式
molecule *this_model;

//! 从模型参数初始化
static PyObject* Model_init(PyObject* self, PyObject* args) {
    para p = default_para;
    if (!PyArg_ParseTuple(args, "iiiiiii", &p.m, &p.n, &p.repeat, 
        &p.direction, &p.glide, &p.climb, &p.bn))
        return NULL;
    int argc = 0; char argv0[] = ""; char *argv = argv0; // 没有参数，待加入
    this_model = new model(p, argc, &argv);
    return Py_None;
}

//! 从文件初始化
static PyObject* Model_restart(PyObject* self, PyObject* args) {
    char fname[64];
    if (!PyArg_ParseTuple(args, "s", &fname))
        return NULL;
    int argc = 0; char argv0[] = ""; char *argv = argv0; // 没有参数，待加入
    this_model = new molecule(std::string(fname), energy_func::node_energy_func,
        energy_func::bond_energy_func, argc, &argv);
    return Py_None;
}

//! 模型结束
static PyObject* Model_delete(PyObject* self, PyObject* args) {
    delete this_model;
    return Py_None;
}

//! 模型更新
static PyObject* Model_update(PyObject* self, PyObject* args) {
    this_model->update();
    return Py_None;
}

//! 模型输出
static PyObject* Model_dump(PyObject* self, PyObject* args) {
    char fname[64];
    nano::dump_t i; 
    if (!PyArg_ParseTuple(args, "si", &fname, &i))
        return NULL;
    this_model->dump(std::string(fname), i);
    return Py_None;
}

//! 总能量输出
static PyObject* Model_energy(PyObject* self, PyObject* args) {
    double up, down;
    if (!PyArg_ParseTuple(args, "dd", &up, &down))
        return NULL;
    double energy = this_model->total_energy(nano::vector(INFINITY, INFINITY,
        down), nano::vector(INFINITY, INFINITY, up));
    return Py_BuildValue("f", &energy);
}

//! 储存到文件
static PyObject* Model_store(PyObject* self, PyObject* args) {
    char fname[64];
    if (!PyArg_ParseTuple(args, "s", &fname))
        return NULL;
    this_model->store(std::string(fname));
    return Py_None;
}

//! 设置参数
static PyObject* Model_set(PyObject* self, PyObject* args) {
    double rest_len, k, tau, precition, step, mass, damp, tempr;  
    if (!PyArg_ParseTuple(args, "dddddddd", &rest_len, &k, &tau, &precition,
        &step, &mass, &damp, &tempr))
        return NULL;
    this_model->set_paras(rest_len, 0); this_model->set_paras(k, 1);
    this_model->set_paras(tau, 2); this_model->set_precision(precition);
    this_model->set_step(step); this_model->set_mass(mass);
    this_model->set_damp(damp); this_model->set_tempr(tempr);
    return Py_None;
}


static PyMethodDef nano_methods[] = {
    {"init", Model_init, METH_VARARGS, "初始化"},
    {"restart", Model_restart, METH_VARARGS, "重启"},
    {"delete", Model_delete, METH_VARARGS, "结束"},
    {"update", Model_update, METH_VARARGS, "更新"},
    {"dump", Model_dump, METH_VARARGS, "输出"},
    {"energy", Model_energy, METH_VARARGS, "能量"},
    {"store", Model_store, METH_VARARGS, "储存"},
    {"set", Model_set, METH_VARARGS, "设置参数"},
    { NULL, NULL}
};
 
static struct PyModuleDef nano_module = {
    PyModuleDef_HEAD_INIT,
    "nano", "", -1,
    nano_methods
};

// 函数名必须与编译出的文件同名
PyMODINIT_FUNC PyInit_nano(void) {
    return PyModule_Create(&nano_module);
}