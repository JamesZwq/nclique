//
// Created by 张文谦 on 25-3-12.
//

#ifndef STATICVECTORE_H
#define STATICVECTORE_H
#include <Global.h>
#include <tbb/spin_mutex.h>

template<typename T>
    class StaticVector {
    public:
        typedef T* iterator;
        typedef const T* const_iterator;

        daf::CliqueSize c_size;
        T *data;

        StaticVector() : c_size(0), data(new T[MAX_CSIZE]) {}

        // copy operator
        StaticVector(const StaticVector &other) : c_size(other.c_size), data(new T[MAX_CSIZE]) {
            memcpy(data, other.data, c_size * sizeof(T));
        }

        StaticVector &operator=(const StaticVector &other) {
            if (this == &other) {
                return *this;
            }
            c_size = other.c_size;
            memcpy(data, other.data, c_size * sizeof(T));
            return *this;
        }

        StaticVector(StaticVector &&other) noexcept : c_size(other.c_size), data(other.data) {
            other.data = nullptr;
            other.c_size = 0;
        }

        operator T *() {
            return data;
        }

        void push_back(const T &value) {
            data[c_size++] = value;
        }

        void pop_back() {
            c_size--;
        }

        bool empty() const {
            return c_size == 0;
        }

        T &operator[](daf::CliqueSize index) {
#if DEBUG
            if (index >= MAX_CSIZE) {
                throw std::out_of_range("Index out of range.");
            }
#endif
            return data[index];
        }

        const T &operator[](daf::CliqueSize index) const {
#if DEBUG
            if (index >= MAX_CSIZE) {
                throw std::out_of_range("Index out of range.");
            }
#endif
            return data[index];
        }

        [[nodiscard]] daf::CliqueSize size() const {
            return c_size;
        }

        // 添加迭代器接口，支持范围 for 循环
        iterator begin() {
            return data;
        }

        iterator end() {
            return data + c_size;
        }

        [[nodiscard]] const_iterator begin() const {
            return data;
        }

        [[nodiscard]] const_iterator end() const {
            return data + c_size;
        }

        [[nodiscard]] const_iterator cbegin() const {
            return data;
        }

        [[nodiscard]] const_iterator cend() const {
            return data + c_size;
        }

        ~StaticVector() {
            delete[] data;
        }
    };


template<typename T>
class MutexStaticVector {
public:
    T *data;
    tbb::spin_mutex *mutexes;

    MutexStaticVector(daf::Size size) : data(new T[size]), mutexes(new tbb::spin_mutex[size]) {}

    void add(daf::Size index, T value) {
        tbb::spin_mutex::scoped_lock lock(mutexes[index]);
        data[index] += value;
    }

    void set(daf::Size index, T value) {
        tbb::spin_mutex::scoped_lock lock(mutexes[index]);
        data[index] = value;
    }

    bool setIfBigger(daf::Size index, T value) {
        tbb::spin_mutex::scoped_lock lock(mutexes[index]);
        if (data[index] < value) {
            data[index] = value;
            return true;
        }
        return false;
    }


};

#endif //STATICVECTORE_H
