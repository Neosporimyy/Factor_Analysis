using System; // Для запуска программы необходимо установить библиотеку MathNet.Numerics
using System.Linq;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

class Program
{
    static void Main()
    {
        // Определение корреляционной матрицы 
        var R = DenseMatrix.OfArray(new double[,]
        {
            {1, 0.92173118157697, 0.599473857228454, 0.821632279232424, -0.200221366165747, -0.461286270601692 },
            {0.92173118157697 , 1,0.562785618658079 ,0.790463543139731 , -0.247144720769157, -0.431121836112182 },
            {0.599473857228454, 0.562785618658079, 1, 0.768403746135087, -0.242237750000099, -0.397312721493995 },
            {0.821632279232424, 0.790463543139731, 0.768403746135087, 1, -0.225948988044967, -0.509849948085866 },
            {-0.200221366165747, -0.247144720769157, -0.242237750000099, -0.225948988044967, 1, 0.722518289024364 },
            {-0.461286270601692, -0.431121836112182, -0.397312721493995, -0.509849948085866,0.722518289024364, 1 }
        });

        // Вычисление собственных значений и собственных векторов 
        var evd = R.Evd();
        List<string> priznak = new List<string> {"Cu","Zn", "V", "Сорг", "Zr", "Ti"};
        List<object> spisok = new List<object>();

        // Собственные значения 
        Console.WriteLine("Собственные значения:");
        foreach (var value in evd.EigenValues)
        {
            spisok.Insert(0, value);
        }
        for (int i = 0; i < spisok.Count; i++)
        {
            Console.WriteLine(priznak[i] + ": " + spisok[i]);
        }

        // Собственные векторы 
        Console.WriteLine("\nСобственные векторы:");
        for (int i = 0; i < evd.EigenVectors.ColumnCount; i++)
        {
            Console.WriteLine($"Собственный вектор {i + 1}:");
            for (int j = 0; j < evd.EigenVectors.RowCount; j++)
            {
                Console.WriteLine($"{j + 1}: {evd.EigenVectors[j, i]}");
            }
            Console.WriteLine();
        }

        // Нормализация собственных векторов
        var normalizedEigenvectors = new DenseMatrix(evd.EigenVectors.RowCount, evd.EigenVectors.ColumnCount);
        List<object> Weight_fact = new List<object>();
        double weight = 0;
        for (int i = 0; i < evd.EigenVectors.ColumnCount; i++)
        {
            var eigenvector = evd.EigenVectors.Column(i);
            var norm = eigenvector.L2Norm();
            normalizedEigenvectors.SetColumn(i, eigenvector / norm);
        }
        int count = 1;
        double eigenvalues = 0;
        int n = 6;
        double summ = 0;
        double[] Eigen_values = new double[n];
        Console.WriteLine("\nФакторные нагрузки:");
        for (int i = evd.EigenVectors.ColumnCount - 1; i >= 0; i--)
        {
            double eigenvalue = evd.EigenValues[i].Real;
            var factorLoadings = normalizedEigenvectors.Column(i) * Math.Sqrt(eigenvalue);
            Console.WriteLine($"Факторные нагрузки для фактора {count}:");
            for (int j = 0; j < factorLoadings.Count; j++)
            {
                Console.WriteLine($"{priznak[j]}: {factorLoadings[j]}");
                weight += factorLoadings[j] * factorLoadings[j];
                summ += factorLoadings[j] * factorLoadings[j];
                eigenvalues += factorLoadings[j];
            }
            Eigen_values[i] = summ;
            summ = 0;
            Console.WriteLine(eigenvalues);
            Console.WriteLine();
            count++;
        }
        Console.WriteLine("Веса факторных нагрузок %");
        int index = 1;
        for(int i = n- 1; i >= 0; i--)
        {
            Console.Write($"F{index}: {Math.Round(Eigen_values[i]/weight*100 , 2)}\n");
            index++;
        }
    }
}
