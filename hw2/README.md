# Домашнее задание №2

- [X] **1. Реализовать классическое перемножение матриц и умножение матрицы на вектор на C/C++ (25%).**

Классическое перемножение реализовано в классе `Matrix`.

- [X] **2. Разбейте на модули, скомпилируйте со статической линковкой, подготовьте Makefile, проверьте флаги -g, -O3 (25%).**

Сборка:
```bash
cd matmul
make
```

Сборка с флагом оптимизации `FLAG`:
```bash
make OPTIMIZE=<FLAG>
```

Выполнение для матриц размера `N`:
```bash
make run N=<N>
```

Очистка:
```bash
make clean
```

- [X] **3. Измерьте времена исполнения для размеров матриц 500, 512, 1000, 1024, 2000, 2048 (25%).**
  * Проведите сравнение с виртуальной машиной, докером (опционально).

Запуск сравнений для классического алгоритма и алгоритма Штрассена с флагами оптимизации `-g` и `-O3` выполняется командой:
```bash
cd matmul
bash timings.sh
```

- [X] **4. Базовые bash-скрипты (25\%).**
```bash
cd bash_scripts
bash 01_for_loop.sh
bash 02_init_print_array.sh
bash 03_calc_float.sh
bash 04_file_exists.sh
```
- [X] **5. Бонус: LINPACK (+20%).**

Запуск бенчмарка выполняется командой:
```bash
cd linpack
bash run.sh
```
Результаты выполнения на локальном компьютере:
```
Sample data file lininput_xeon64.

Current date/time: Fri Oct  8 20:24:34 2021

CPU frequency:    3.697 GHz
Number of CPUs: 1
Number of cores: 6
Number of threads: 6

Parameters are set to:

Number of tests: 5
Number of equations to solve (problem size) : 1000  2000  5000  10000 20000
Leading dimension of array                  : 1000  2000  5008  10000 20000
Number of trials to run                     : 4     2     2     2     1    
Data alignment value (in Kbytes)            : 4     4     4     4     4    

Maximum memory requested that can be used=3200404096, at the size=20000

=================== Timing linear equation system solver ===================

Size   LDA    Align. Time(s)    GFlops   Residual     Residual(norm) Check
1000   1000   4      0.007      99.4748  1.160572e-12 3.469815e-02   pass
1000   1000   4      0.005      121.9759 1.160572e-12 3.469815e-02   pass
1000   1000   4      0.006      118.4922 1.160572e-12 3.469815e-02   pass
1000   1000   4      0.006      111.9508 1.160572e-12 3.469815e-02   pass
2000   2000   4      0.039      138.2230 3.396561e-12 2.702137e-02   pass
2000   2000   4      0.032      165.8536 3.396561e-12 2.702137e-02   pass
5000   5008   4      0.484      172.2796 2.230932e-11 2.943595e-02   pass
5000   5008   4      0.499      167.0092 2.230932e-11 2.943595e-02   pass
10000  10000  4      2.847      234.2133 1.014424e-10 3.428008e-02   pass
10000  10000  4      2.538      262.7777 1.014424e-10 3.428008e-02   pass
20000  20000  4      19.031     280.2861 3.811837e-10 3.251820e-02   pass

Performance Summary (GFlops)

Size   LDA    Align.  Average  Maximal
1000   1000   4       112.9734 121.9759
2000   2000   4       152.0383 165.8536
5000   5008   4       169.6444 172.2796
10000  10000  4       248.4955 262.7777
20000  20000  4       280.2861 280.2861

Residual checks PASSED

End of tests
```

- [X] **6. Супербонус: протестируйте алгоритм Штрассена (+20%).**

Алгоритм Штрассена реализован в функции `strassen_product` (см. `matrix.cpp`). Тестирование производительности производится совместно с классическим алгоритмом (см. пункт 3).