// Copyright 2021 Dobrov Pavel
#include <gtest/gtest.h>
// #include <chrono>
#include <vector>
#include "./m_cannon.h"

// TEST(m_Cannon, multi_1000x1000) {
//     std::vector<double> A = RandomMatrix(1000);
//     std::vector<double> B = RandomMatrix(1000);

//     auto start = std::chrono::high_resolution_clock::now();
//     std::vector<double> C1 = seqMulti(A, B, 1000);
//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double, std::milli> duration = end - start;
//     printf("Seq time: %lf s\n", duration.count() / 1000);
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<double> C2 = parMulti(A, B, 1000);
//     end = std::chrono::high_resolution_clock::now();
//     duration = end - start;
//     printf("Par time: %lf s\n", duration.count() / 1000);

//     ASSERT_EQ(C1, C2);
// }

TEST(m_Cannon, matrix_with_negative_size) {
    ASSERT_ANY_THROW(RandomMatrix(-1));
}

TEST(m_Cannon, Cant_execute_with_negative_size_parMulti) {
    std::vector<double> A(1);
    std::vector<double> B(1);
    ASSERT_ANY_THROW(parMulti(A, B, -1));
}

TEST(m_Cannon, milti_123x123) {
    std::vector<double> A = RandomMatrix(123);
    std::vector<double> B = RandomMatrix(123);
    std::vector<double> C1 = seqMulti(A, B, 123);
    std::vector<double> C2 = parMulti(A, B, 123);

    ASSERT_EQ(C1, C2);
}

TEST(m_Cannon, milti_1x1) {
    std::vector<double> A = RandomMatrix(1);
    std::vector<double> B = RandomMatrix(1);
    std::vector<double> C1 = seqMulti(A, B, 1);
    std::vector<double> C2 = parMulti(A, B, 1);

    ASSERT_EQ(C1, C2);
}

TEST(m_Cannon, milti_101x101) {
    std::vector<double> A = RandomMatrix(101);
    std::vector<double> B = RandomMatrix(101);
    std::vector<double> C1 = seqMulti(A, B, 101);
    std::vector<double> C2 = parMulti(A, B, 101);

    ASSERT_EQ(C1, C2);
}


int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}
