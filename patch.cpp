// c++ version of patching program

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif
using std::cerr;
using std::clock;
using std::clock_t;
using std::cout;
using std::endl;
using std::ifstream;
using std::map;
using std::ofstream;
using std::setfill;
using std::setprecision;
using std::setw;
using std::string;
using std::vector;

// #define nhist 2
// #define ncompin 2

// Gaussian elimination for Ax = b
// Solves in-place, modifying A and b. Returns false if singular.
bool gauss_solve(std::vector<std::vector<double>> &A, std::vector<double> &b,
                 std::vector<double> &x) {
  int n = A.size();
  for (int i = 0; i < n; i++) {
    double maxEl = fabs(A[i][i]);
    int maxRow = i;
    for (int k = i + 1; k < n; k++) {
      if (fabs(A[k][i]) > maxEl) {
        maxEl = fabs(A[k][i]);
        maxRow = k;
      }
    }
    if (maxEl < 1e-20)
      return false;

    for (int k = i; k < n; k++)
      std::swap(A[maxRow][k], A[i][k]);
    std::swap(b[maxRow], b[i]);

    for (int k = i + 1; k < n; k++) {
      double c = -A[k][i] / A[i][i];
      for (int j = i; j < n; j++) {
        if (i == j)
          A[k][j] = 0;
        else
          A[k][j] += c * A[i][j];
      }
      b[k] += c * b[i];
    }
  }
  x.assign(n, 0.0);
  for (int i = n - 1; i >= 0; i--) {
    x[i] = b[i];
    for (int k = i + 1; k < n; k++)
      x[i] -= A[i][k] * x[k];
    x[i] = x[i] / A[i][i];
  }
  return true;
}

int main(int argc, char *argv[]) {
  // Parse command line arguments for OpenMP threads
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '+' && argv[i][1] == 'p') {
      int num_threads = atoi(&argv[i][2]);
      if (num_threads > 0) {
#ifdef _OPENMP
        omp_set_num_threads(num_threads);
        cout << "Setting OpenMP threads to " << num_threads << endl;
#else
        cout << "OpenMP not enabled during compilation. Ignoring thread setting." << endl;
#endif
      }
    }
  }

  // read histogram labels
  clock_t start_time = std::clock();
  string suffix = "a";
  int ncompin;
  int nhist;

  ifstream fin;
  fin.open("input_hs.dat");
  if (!fin.is_open()) {
    cerr << "Error: Could not open input_hs.dat" << endl;
    exit(1);
  }

  vector<int> histid;

  fin >> ncompin;
  fin >> suffix;
  fin >> nhist;

  int temp;

  for (int i = 0; i < nhist; ++i) {
    fin >> temp;
    histid.push_back(temp);
  }

  fin.close();

  // build filenames and read histograms

  double *t = new double[nhist];
  double *beta = new double[nhist]; // inverse temperature
  int *ncomp = new int[nhist];

  double **mu = new double *[nhist];
  for (int i = 0; i < nhist; ++i) {
    mu[i] = new double[ncompin];
  }

  double mu1;
  double mu2;
  double *weight = new double[nhist];
  double *oldweight = new double[nhist];
  double *lx = new double[nhist];
  double *ly = new double[nhist];
  double *lz = new double[nhist];

  int **minp = new int *[nhist];
  for (int i = 0; i < nhist; ++i) {
    minp[i] = new int[ncompin];
  }

  int **maxp = new int *[nhist];
  for (int i = 0; i < nhist; ++i) {
    maxp[i] = new int[ncompin];
  }

  int *maxp_overall = new int[ncompin];
  for (int j = 0; j < ncompin; ++j) {
    maxp_overall[j] = -1e5;
  }
  vector<string> histname;
  double xbox;
  double ybox;
  double zbox;
  vector<vector<vector<int>>> n(ncompin);
  vector<vector<double>> e;
  vector<vector<int>> nt_temp(ncompin);
  vector<double> et;
  vector<int> nentry;
  double etemp;
  double crap;
  string s;
  double a;
  double b;
  int ncompin_old;
  int nhist_old;
  bool oldrun;

  // see if we already have weights to use as an initial guess
  ifstream fin2;
  fin2.open("weights.dat");
  if (fin2) {
    oldrun = true;
    cout << "Reading weights from: weights.dat" << endl;
    cout << setfill('=') << setw(60) << "=" << setfill(' ') << endl;

    int i = 0;
    string temp;
    fin2 >> ncompin_old;
    fin2 >> nhist_old;
    fin2 >> suffix;
    for (int i = 0; i < nhist_old; ++i) {
      fin2 >> histid[i] >> temp >> oldweight[i] >> t[i];
      cout << setw(4) << histid[i] << std::scientific << setprecision(6)
           << setw(14) << oldweight[i] << std::fixed << setprecision(2)
           << setw(8) << t[i];
      for (int c = 0; c < ncompin; ++c) {
        fin2 >> mu[i][c];
        cout << setw(10) << mu[i][c];
      }
      cout << endl;
    }
    fin2.clear();
    fin2.close();
    if (nhist > nhist_old) {
      for (int i = nhist_old; i < nhist; ++i) {
        oldweight[i] = 100000000.0;
      }
    }

  } else if (!fin2) {
    oldrun = false;
    for (int i = 0; i < nhist; ++i) {
      oldweight[i] = 100000000.0;
    }
  }
  cout << setfill('=') << setw(60) << "=" << setfill(' ') << endl;

  cout << "Patching histograms" << endl;
  cout << ncompin << " component system" << endl;
  cout << "Reading " << nhist << " histograms" << endl;
  cout << endl;
  // cout << setfill('=')<<setw(80) << "=" <<setfill(' ')<<endl;
  cout << setw(8) << "Histogram" << setw(8) << "Nentry" << setw(8) << "Temp";
  for (int c = 0; c < ncompin; ++c) {
    cout << setw(8) << "mu" << c + 1;
  }
  for (int c = 0; c < ncompin; ++c) {
    cout << setw(8) << "n" << c + 1 << "-min" << setw(8) << "n" << c + 1 << "-max";
  }
  cout << endl;
  cout << setfill('=') << setw(80) << "=" << setfill(' ') << endl;
  cout << std::fixed;
  for (int i = 0; i < nhist; ++i) {
    // build file names
    string filein;
    for (int j = 0; j < ncompin; ++j) {
      minp[i][j] = 1e5;
      maxp[i][j] = -1e5;
    }
    filein = "his" + std::to_string(histid[i]) + suffix + ".dat";
    histname.push_back(filein);
    fin.open(filein.c_str()); // You can pass a string
    if (!fin.is_open()) {
      cerr << "Error: Could not open histogram file " << filein << endl;
      exit(1);
    }
    // read header information
    fin >> t[i] >> ncomp[i];
    for (int c = 0; c < ncompin; ++c) {
      fin >> mu[i][c];
    }
    fin >> lx[i] >> ly[i] >> lz[i];

    beta[i] = 1.0 / t[i]; // Computing the reciprocal value for the temp

    // error checking
    if (ncompin != ncomp[i]) {
      cerr << "Mismatch between number of components in input_hs.dat and "
              "histogram file "
           << histname[i] << endl;
      cerr << "Components in input_hs.dat = " << ncompin << endl;
      cerr << "Components in histogram = " << ncomp[i] << endl;
      exit(1);
    }
    // set initial box size.  Check all other histograms to make sure they were
    // run for the same system size.
    if (i == 0) {
      xbox = lx[i];
      ybox = ly[i];
      zbox = lz[i];
    } else if (i != 0) {
      if (xbox != lx[i] || ybox != ly[i] || zbox != lz[i]) {
        cerr << "System size in histogram " << filein
             << " is inconsistent with other files" << endl;
        exit(1);
      }
    }
    int j = 0;
    vector<int> ntemp(ncompin);
    while (true) {
      for (int c = 0; c < ncompin; ++c) {
        fin >> ntemp[c];
      }
      fin >> etemp;
      if (fin.eof() || fin.fail()) break;

      for (int c = 0; c < ncompin; ++c) {
        minp[i][c] = std::min(minp[i][c], ntemp[c]);
        maxp[i][c] = std::max(maxp[i][c], ntemp[c]);
        maxp_overall[c] = std::max(maxp[i][c], maxp_overall[c]);
        nt_temp[c].push_back(ntemp[c]);
      }
      et.push_back(etemp);
      ++j;
    }
    // move data from temporary variables into vectors.  First index is
    // component,
    //  second index is file name, third is the data label
    for (int c = 0; c < ncompin; ++c) {
      n[c].push_back(nt_temp[c]);
      nt_temp[c].clear();
    }
    e.push_back(et);
    nentry.push_back(n[0][i].size());
    cout << setprecision(2) << setw(8) << histname[i] << setw(8) << nentry[i]
         << setw(8) << t[i];
    for (int c = 0; c < ncompin; ++c) {
      cout << setw(10) << mu[i][c];
    }
    for (int c = 0; c < ncompin; ++c) {
      cout << setw(8) << minp[i][c] << setw(8) << maxp[i][c];
    }
    cout << endl;

    // clear temporary arrays, otherwise "bad things" happen
    et.clear();

    fin.clear();
    fin.close();
  }
  cout << "Histogram data read successfully" << endl;

  int iter = 0;
  double maxd = 100;
  double tol = 5e-5;
  double val;
  int N_total = 0;

  for (int i = 0; i < nhist; ++i) {
    N_total += nentry[i];
  }

  // Precalculate A_j(data) = -beta_j * E + sum_c (beta_j * mu_{j,c} * N_c)
  // This replaces the enormous 'prob' array that caused massive over-allocation
  double *A = new double[N_total * nhist];
  int data_idx = 0;
  for (int ifile = 0; ifile < nhist; ++ifile) {
    for (int i = 0; i < nentry[ifile]; ++i) {
      for (int jfile = 0; jfile < nhist; ++jfile) {
        double a_val = -beta[jfile] * e[ifile][i];
        for (int icomp = 0; icomp < ncompin; ++icomp) {
          a_val += beta[jfile] * mu[jfile][icomp] * n[icomp][ifile][i];
        }
        A[data_idx * nhist + jfile] = a_val;
      }
      data_idx++;
    }
  }

  // Initialize weights using Mean-Energy Approximation for histograms without a
  // starting guess
  int d_start = 0;
  vector<double> log_w(nhist, 0.0);
  bool has_new_weights = false;

  for (int jfile = 0; jfile < nhist; ++jfile) {
    if (oldweight[jfile] == 100000000.0) {
      double avg_A = 0.0;
      for (int i = 0; i < nentry[jfile]; ++i) {
        avg_A += A[(d_start + i) * nhist + jfile];
      }
      if (nentry[jfile] > 0)
        avg_A /= nentry[jfile];
      log_w[jfile] = avg_A;
      has_new_weights = true;
    }
    d_start += nentry[jfile];
  }

  if (has_new_weights) {
    double max_log = -1e9;
    for (int jfile = 0; jfile < nhist; ++jfile) {
      if (oldweight[jfile] == 100000000.0 && log_w[jfile] > max_log) {
        max_log = log_w[jfile];
      }
    }
    for (int jfile = 0; jfile < nhist; ++jfile) {
      if (oldweight[jfile] == 100000000.0) {
        // exp(avg_A - max) ensures it is bounded between 0 and 1
        oldweight[jfile] = exp(log_w[jfile] - max_log);
      }
    }
  }

  double *log_N_over_f = new double[nhist];
  double *denom_log = new double[N_total];

  // DIIS parameters and buffers
  const int diis_m = 5;
  std::vector<std::vector<double>> diis_errs;
  std::vector<std::vector<double>> diis_weights;

  while (maxd > tol && iter < 1e4) {
    ++iter;
    maxd = 0.0;

    // 1. Precalculate loop invariants for this iteration
    for (int jfile = 0; jfile < nhist; ++jfile) {
      log_N_over_f[jfile] = log(nentry[jfile] / oldweight[jfile]);
    }

// 2. Compute denominator log-sum-exp for each data point
#pragma omp parallel for
    for (int d = 0; d < N_total; ++d) {
      double ylog = -1e9;
      for (int jfile = 0; jfile < nhist; ++jfile) {
        double tmp1 = A[d * nhist + jfile] + log_N_over_f[jfile];
        ylog = std::max(ylog, tmp1) + log(1.0 + exp(-fabs(ylog - tmp1)));
      }
      denom_log[d] = ylog;
    }

    std::vector<double> current_weight(nhist, 0.0);

    int nthreads = 1;
#pragma omp parallel
    {
#pragma omp single
#ifdef _OPENMP
      nthreads = omp_get_num_threads();
#else
      nthreads = 1;
#endif
    }

    std::vector<std::vector<double>> thread_weights(
        nthreads, std::vector<double>(nhist, 0.0));

#pragma omp parallel
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      // Use a purely local array to prevent False Sharing (cache line
      // contention) between threads
      std::vector<double> local_weights(nhist, 0.0);

#pragma omp for schedule(static)
      for (int d = 0; d < N_total; ++d) {
        for (int kfile = 0; kfile < nhist; ++kfile) {
          local_weights[kfile] += exp(A[d * nhist + kfile] - denom_log[d]);
        }
      }

      // Write back to the shared array only once at the very end
      for (int kfile = 0; kfile < nhist; ++kfile) {
        thread_weights[tid][kfile] = local_weights[kfile];
      }
    }

    // Deterministic reduction to prevent floating-point noise from
    // destabilizing DIIS
    for (int t = 0; t < nthreads; ++t) {
      for (int kfile = 0; kfile < nhist; ++kfile) {
        current_weight[kfile] += thread_weights[t][kfile];
      }
    }

    // Calculate errors for DIIS
    std::vector<double> current_err(nhist, 0.0);
    for (int i = 0; i < nhist; ++i) {
      current_err[i] = current_weight[i] - oldweight[i];
    }

    // Update DIIS buffers
    diis_weights.push_back(current_weight);
    diis_errs.push_back(current_err);
    if (diis_weights.size() > diis_m) {
      diis_weights.erase(diis_weights.begin());
      diis_errs.erase(diis_errs.begin());
    }

    bool diis_used = false;
    std::vector<double> extrap_weight(nhist, 0.0);

    // Attempt DIIS extrapolation
    int n_diis = diis_errs.size();
    if (n_diis >= 3) {
      std::vector<std::vector<double>> B(n_diis + 1,
                                         std::vector<double>(n_diis + 1, 0.0));
      std::vector<double> b_vec(n_diis + 1, 0.0);
      std::vector<double> c_vec(n_diis + 1, 0.0);

      for (int i = 0; i < n_diis; ++i) {
        for (int j = 0; j <= i; ++j) {
          double dot = 0.0;
          for (int k = 0; k < nhist; ++k)
            dot += diis_errs[i][k] * diis_errs[j][k];
          B[i][j] = B[j][i] = dot;
        }
        B[i][n_diis] = B[n_diis][i] = -1.0;
      }
      B[n_diis][n_diis] = 0.0;
      b_vec[n_diis] = -1.0;

      if (gauss_solve(B, b_vec, c_vec)) {
        for (int i = 0; i < n_diis; ++i) {
          for (int k = 0; k < nhist; ++k) {
            extrap_weight[k] += c_vec[i] * diis_weights[i][k];
          }
        }
        bool valid = true;
        for (int k = 0; k < nhist; ++k) {
          if (extrap_weight[k] <= 0.0 || std::isnan(extrap_weight[k]))
            valid = false;
        }
        if (valid)
          diis_used = true;
      }
    }

    // Determine maxd using the non-extrapolated fixed point step for true
    // convergence
    for (int kfile = 0; kfile < nhist; ++kfile) {
      val = fabs(current_weight[kfile] / oldweight[kfile] - 1.0);
      if (val > maxd)
        maxd = val;
    }

    // Assign the new weight, using extrapolation if successful
    for (int kfile = 0; kfile < nhist; ++kfile) {
      weight[kfile] = diis_used ? extrap_weight[kfile] : current_weight[kfile];
    }

    // Clear DIIS if diverged numerically
    if (diis_used && maxd > 10.0) {
      diis_weights.clear();
      diis_errs.clear();
      for (int kfile = 0; kfile < nhist; ++kfile)
        weight[kfile] = current_weight[kfile];
      maxd = 100.0; // Force another iter
    }

    // normalize weights
    for (int i = 1; i < nhist; ++i) {
      weight[i] = weight[i] / weight[0];
      oldweight[i] = weight[i];
    }
    weight[0] = 1.0;
    oldweight[0] = weight[0];

    // Conditional output: write every 10 iterations or at convergence
    if (iter % 10 == 0 || maxd <= tol) {
      ofstream fout;
      fout.open("weights.dat");
      if (!fout.is_open()) {
        cerr << "Error: Could not open weights.dat for writing." << endl;
        exit(1);
      }
      fout << ncompin << endl;
      fout << nhist << endl;
      fout << suffix << endl;
      for (int ifile = 0; ifile < nhist; ++ifile) {
        fout << setw(4) << histid[ifile] << setw(10) << nentry[ifile]
             << setw(20) << std::scientific << setprecision(8) << weight[ifile]
             << setw(10) << std::fixed << setprecision(2) << t[ifile];
        for (int c = 0; c < ncompin; ++c) {
          fout << setw(10) << mu[ifile][c];
        }
        fout << endl;
      }
      fout << "Total iterations = " << iter << endl;
      fout << "Convergence = " << setprecision(10) << maxd << endl;
      fout << "Run time " << ((clock() - start_time) / (double)CLOCKS_PER_SEC)
           << " (s)" << endl;
      fout.close();
    }
  }

  delete[] A;
  delete[] log_N_over_f;
  delete[] denom_log;

  return 0;
}
