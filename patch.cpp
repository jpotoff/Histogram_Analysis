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

  for (int i = 0; i < nhist; ++i) {
    for (int j = 0; j < i; ++j) {
      if (histid[i] == histid[j]) {
        cerr << "Error: histogram " << histid[i]
             << " is listed more than once in input_hs.dat" << endl;
        exit(1);
      }
    }
  }

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

  // see if we already have weights to use as an initial guess.  The
  // histograms in weights.dat must be the first entries of input_hs.dat, in
  // the same order; new histograms may only be appended to the end.
  ifstream fin2;
  fin2.open("weights.dat");
  if (fin2) {
    oldrun = true;
    cout << "Reading weights from: weights.dat" << endl;
    cout << setfill('=') << setw(60) << "=" << setfill(' ') << endl;

    string suffix_old;
    fin2 >> ncompin_old;
    fin2 >> nhist_old;
    fin2 >> suffix_old;
    if (!fin2) {
      cerr << "Error: could not read the header of weights.dat" << endl;
      exit(1);
    }
    if (ncompin_old != ncompin) {
      cerr << "Error: weights.dat is for a " << ncompin_old
           << " component system, but input_hs.dat specifies " << ncompin
           << " components" << endl;
      exit(1);
    }
    if (suffix_old != suffix) {
      cerr << "Error: weights.dat was generated for suffix '" << suffix_old
           << "', but input_hs.dat specifies suffix '" << suffix << "'"
           << endl;
      exit(1);
    }
    if (nhist_old > nhist) {
      cerr << "Error: weights.dat contains " << nhist_old
           << " histograms, but input_hs.dat lists only " << nhist << endl;
      cerr << "Histograms may only be added to the end of input_hs.dat. "
              "Delete weights.dat to start from scratch."
           << endl;
      exit(1);
    }
    string temp;
    for (int i = 0; i < nhist_old; ++i) {
      int id_old;
      fin2 >> id_old >> temp >> oldweight[i] >> t[i];
      for (int c = 0; c < ncompin; ++c) {
        fin2 >> mu[i][c];
      }
      if (!fin2) {
        cerr << "Error: could not read entry " << i + 1 << " of weights.dat"
             << endl;
        exit(1);
      }
      if (id_old != histid[i]) {
        cerr << "Error: histogram order in input_hs.dat does not match "
                "weights.dat"
             << endl;
        cerr << "Entry " << i + 1 << " is histogram " << histid[i]
             << " in input_hs.dat but histogram " << id_old
             << " in weights.dat" << endl;
        cerr << "The histograms in weights.dat must be listed first in "
                "input_hs.dat, in the same order;"
             << endl;
        cerr << "new histograms may only be added to the end. Fix "
                "input_hs.dat, or delete weights.dat to start from scratch."
             << endl;
        exit(1);
      }
      cout << setw(4) << histid[i] << std::scientific << setprecision(6)
           << setw(14) << oldweight[i] << std::fixed << setprecision(2)
           << setw(8) << t[i];
      for (int c = 0; c < ncompin; ++c) {
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

  for (int jfile = 0; jfile < nhist; ++jfile) {
    if (oldweight[jfile] == 100000000.0) {
      double avg_A = 0.0;
      for (int i = 0; i < nentry[jfile]; ++i) {
        avg_A += A[(d_start + i) * nhist + jfile];
      }
      if (nentry[jfile] > 0)
        avg_A /= nentry[jfile];
      log_w[jfile] = avg_A;
          }
    d_start += nentry[jfile];
  }

  // Histograms without a starting weight get ln(weight) = <A>_j, shifted by
  // the largest such value.  This is kept in log form: for large systems
  // ln(weight) spans thousands and exp() would overflow or underflow.
  double max_log = -1e300;
  for (int jfile = 0; jfile < nhist; ++jfile) {
    if (oldweight[jfile] == 100000000.0 && log_w[jfile] > max_log) {
      max_log = log_w[jfile];
    }
  }

  // Solve the self-consistent histogram equations by Newton-Raphson in log
  // space.  With g_j = ln(weight_j), the weights are the minimizer of the
  // convex function (Shirts & Chodera, J. Chem. Phys. 129, 124105 (2008))
  //   Phi(g) = sum_d log sum_j N_j exp(A_dj - g_j) + sum_j N_j g_j
  // whose gradient is grad_k = N_k - W_k, with W_k / N_k = (fixed-point
  // update) / (current weight).  grad = 0 is the same solution as the old
  // iteration, and maxd = max_k |W_k / N_k - 1| is the same convergence test.
  // Working with ln(weight) avoids the many decades spanned by the weights.
  // g_0 = 0 fixes the arbitrary normalization.
  vector<double> g(nhist), log_N(nhist);
  for (int j = 0; j < nhist; ++j) {
    g[j] = (oldweight[j] == 100000000.0) ? log_w[j] - max_log
                                          : log(oldweight[j]);
    log_N[j] = log((double)nentry[j]);
  }
  for (int j = nhist - 1; j >= 0; --j)
    g[j] -= g[0];
  double *denom_log = new double[N_total];

  // Evaluates Phi(gg) and fills denom_log[d] = log sum_j N_j exp(A_dj - gg_j)
  auto objective = [&](const vector<double> &gg) {
    double phi = 0.0;
#pragma omp parallel for reduction(+ : phi)
    for (int d = 0; d < N_total; ++d) {
      double ymax = -1e300;
      for (int j = 0; j < nhist; ++j)
        ymax = std::max(ymax, A[d * nhist + j] + log_N[j] - gg[j]);
      double sum = 0.0;
      for (int j = 0; j < nhist; ++j)
        sum += exp(A[d * nhist + j] + log_N[j] - gg[j] - ymax);
      denom_log[d] = ymax + log(sum);
      phi += denom_log[d];
    }
    for (int j = 0; j < nhist; ++j)
      phi += nentry[j] * gg[j];
    return phi;
  };

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

  const int max_iter = 1000;
  const double max_step = 1000.0; // largest allowed change in any ln(weight)
  const double step_tol = 1e-6;   // converged ln(weight) changes less than this
  double phi = objective(g);
  vector<double> W(nhist), grad(nhist), H(nhist * nhist);

  while (true) {
    ++iter;

    // W_k = sum_d p_dk and Hessian H_kl = sum_d (p_dk delta_kl - p_dk p_dl),
    // with p_dk = N_k exp(A_dk - g_k) / denom_d.  Per-thread accumulation
    // followed by a deterministic reduction.
    vector<vector<double>> thread_W(nthreads, vector<double>(nhist, 0.0));
    vector<vector<double>> thread_H(nthreads,
                                    vector<double>(nhist * nhist, 0.0));
#pragma omp parallel
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      vector<double> lW(nhist, 0.0), lH(nhist * nhist, 0.0), p(nhist);
#pragma omp for schedule(static)
      for (int d = 0; d < N_total; ++d) {
        for (int k = 0; k < nhist; ++k) {
          p[k] = exp(A[d * nhist + k] + log_N[k] - g[k] - denom_log[d]);
          lW[k] += p[k];
        }
        // Diagonal written as p_k * sum_{j!=k} p_j rather than p_k - p_k^2,
        // which cancels catastrophically when one histogram has p_k ~ 1
        for (int k = 0; k < nhist; ++k)
          for (int l = 0; l < nhist; ++l)
            if (l != k) {
              lH[k * nhist + l] -= p[k] * p[l];
              lH[k * nhist + k] += p[k] * p[l];
            }
      }
      thread_W[tid] = lW;
      thread_H[tid] = lH;
    }
    std::fill(W.begin(), W.end(), 0.0);
    std::fill(H.begin(), H.end(), 0.0);
    for (int t = 0; t < nthreads; ++t) {
      for (int k = 0; k < nhist; ++k)
        W[k] += thread_W[t][k];
      for (int kl = 0; kl < nhist * nhist; ++kl)
        H[kl] += thread_H[t][kl];
    }
    maxd = 0.0;
    for (int k = 0; k < nhist; ++k) {
      grad[k] = nentry[k] - W[k];
      maxd = std::max(maxd, fabs(W[k] / nentry[k] - 1.0));
    }
    if (!std::isfinite(maxd)) {
      cerr << "Error: non-finite weights at iteration " << iter << endl;
      exit(1);
    }

    // Newton step on the reduced system (g_0 held fixed)
    int nr = nhist - 1;
    vector<vector<double>> Hr(nr, vector<double>(nr));
    vector<double> rhs(nr), step_r;
    for (int k = 0; k < nr; ++k) {
      rhs[k] = -grad[k + 1];
      for (int l = 0; l < nr; ++l)
        Hr[k][l] = H[(k + 1) * nhist + (l + 1)];
    }
    vector<double> step(nhist, 0.0);
    bool newton_ok = gauss_solve(Hr, rhs, step_r);
    double max_newton = 0.0;
    if (newton_ok) {
      for (int k = 0; k < nr; ++k) {
        step[k + 1] = step_r[k];
        max_newton = std::max(max_newton, fabs(step_r[k]));
      }
    }

    // A small maxd alone is not enough: a weakly overlapping histogram can
    // still be off in ln(weight), so also require the Newton step to vanish
    bool done = (maxd <= tol && newton_ok && max_newton <= step_tol) ||
                iter >= max_iter;

    cout << "Iteration " << setw(5) << iter << "   convergence = "
         << std::scientific << setprecision(4) << maxd;
    if (newton_ok)
      cout << "   Newton step = " << max_newton;
    cout << std::fixed << endl;

    // Conditional output: write every 10 iterations or at the end
    if (iter % 10 == 0 || done) {
      for (int j = 0; j < nhist; ++j)
        weight[j] = exp(g[j]);
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
    if (done)
      break;

    // Backtracking line search on Phi along step (Armijo condition).
    // Updates g, phi and denom_log and returns true on success.
    vector<double> g_new(nhist);
    auto line_search = [&](vector<double> step) {
      double big = 0.0;
      for (int k = 0; k < nhist; ++k)
        big = std::max(big, fabs(step[k]));
      if (big > max_step)
        for (int k = 0; k < nhist; ++k)
          step[k] *= max_step / big;
      double slope = 0.0;
      for (int k = 0; k < nhist; ++k)
        slope += grad[k] * step[k];
      if (!(slope < 0.0))
        return false; // not a descent direction
      double lambda = 1.0;
      for (int ls = 0; ls < 40; ++ls) {
        for (int k = 0; k < nhist; ++k)
          g_new[k] = g[k] + lambda * step[k];
        double phi_new = objective(g_new);
        if (phi_new <= phi + 1e-4 * lambda * slope + 1e-12 * fabs(phi)) {
          g = g_new;
          phi = phi_new; // denom_log now corresponds to g
          return true;
        }
        lambda *= 0.5;
      }
      objective(g); // restore denom_log for the unchanged g
      return false;
    };

    bool moved = false;
    if (newton_ok)
      moved = line_search(step);
    if (!moved) {
      // Newton failed (nearly singular Hessian far from the solution, where
      // histograms do not yet overlap): take a fixed-point step instead
      for (int k = 1; k < nhist; ++k)
        step[k] = log(std::max(W[k], 1e-300) / nentry[k]) -
                  log(std::max(W[0], 1e-300) / nentry[0]);
      moved = line_search(step);
    }
    if (!moved) {
      cerr << "Error: line search failed at iteration " << iter << endl;
      break;
    }
  }

  if (iter >= max_iter) {
    cerr << "Warning: weights did not converge in " << max_iter
         << " iterations (convergence = " << std::scientific << maxd << ")"
         << endl;
  }

  delete[] A;
  delete[] denom_log;

  return 0;
}
