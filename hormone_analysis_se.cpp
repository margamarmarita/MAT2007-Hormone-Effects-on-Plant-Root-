// Libraries needed for downloading, reading, analysing and saving the data
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <cstdlib>
using namespace std;

// MAT2007 Project - Statistical Analysis of Root-Hormone: Effect of Phytohormones on Plant 
//                   Root Lengths and the Correlations of Change



struct Condition {
    string label;               // Name of the treatment, such as IAA or ABA
    int phenotype_id;           // AraPheno ID used to find the correct dataset
    map<int, double> values;    // Root length for each accession
    string pheno_name;          // Phenotype name read from the downloaded file
};


// EXERCISE 1 - Download the data and load it into the programme

// Separates a CSV line into individual columns - commas inside quotation
// marks are ignored so that they do not split a field by mistake
vector<string> split_csv_line(const string& line) {
    vector<string> fields;
    string current;
    bool in_quotes = false;

    for (char c : line) {
        if (c == '"') {
            in_quotes = !in_quotes;
        } else if (c == ',' && !in_quotes) {
            fields.push_back(current);
            current.clear();
        } else if (c != '\r') {
            current.push_back(c);
        }
    }

    fields.push_back(current);
    return fields;
}

// Checks whether the CSV file already exists - if not, curl is used
// to download it from AraPheno.
bool download_if_missing(int phenotype_id, const string& filename) {
    ifstream test(filename);
    if (test) return true; // No need to download the file again

    string url = "https://arapheno.1001genomes.org/rest/phenotype/" +
                 to_string(phenotype_id) + "/values.csv";
    string command = "curl -L -s -o " + filename + " " + url;

    cout << "Downloading " << url << endl;
    int status = system(command.c_str());

    // Check that the download worked and that the file can be opened
    ifstream check(filename);
    return (status == 0 && check);
}

// Reads the CSV file and stores each measurement alongside its accession ID
//      The column names are used to find the right data, so their order does not matter
bool load_phenotype(const string& filename, map<int, double>& data,
                    string& pheno_name) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: Could not open " << filename << endl;
        return false;
    }

    string line;
    getline(infile, line); // Read the column names from the first line
    vector<string> header = split_csv_line(line);

    int col_id = -1, col_val = -1, col_name = -1;

    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == "accession_id") col_id = (int)i;
        else if (header[i] == "phenotype_value") col_val = (int)i;
        else if (header[i] == "phenotype_name") col_name = (int)i;
    }

    // Stop if the file does not contain the columns needed for the analysis
    if (col_id < 0 || col_val < 0) {
        cerr << "Error: columns 'accession_id' / 'phenotype_value' not found in "
             << filename << ". Header was: " << line << endl;
        return false;
    }

    while (getline(infile, line)) {
        if (line.empty()) continue;

        vector<string> f = split_csv_line(line);
        if (f.size() <= (size_t)max(col_id, col_val)) continue;
        if (f[col_val].empty()) continue; // Skip accessions without a measurement

        try {
            int id = stoi(f[col_id]);
            double v = stod(f[col_val]);
            data[id] = v;

            if (col_name >= 0 && f.size() > (size_t)col_name)
                pheno_name = f[col_name];
        } catch (...) {
            continue; // Ignore rows containing values that cannot be read
        }
    }

    return !data.empty();
}


// EXERCISE 2 - Functions used to calculate the statistics

// Makes a shuffled copy of the measurements so that the original data
//      stays unchanged when randomly dividing the accessions into groups
vector<double> shuffled_copy(const vector<double>& data, mt19937& gen) {
    vector<double> copy = data;
    shuffle(copy.begin(), copy.end(), gen);
    return copy;
}

// Calculates the average value in a dataset
double mean(const vector<double>& data) {
    return accumulate(data.begin(), data.end(), 0.0) / data.size();
}

// Calculates the population standard deviation, using N in the denominator
//      This follows the approach used in the course tutorials
double standard_deviation(const vector<double>& data) {
    double m = mean(data);
    double sum_sq = 0.0;

    for (double x : data)
        sum_sq += (x - m) * (x - m);

    return sqrt(sum_sq / data.size());
}

// Calculates the standard error of the mean using the sample standard
//      deviation and all the available measurements
double standard_error(const vector<double>& data) {
    double m = mean(data);
    double sum_sq = 0.0;

    for (double x : data)
        sum_sq += (x - m) * (x - m);

    double sample_std = sqrt(sum_sq / (data.size() - 1));
    return sample_std / sqrt((double)data.size());
}

// Estimates the standard error of a Pearson correlation coefficient
//      from the correlation value and the number of paired measurements
double pearson_standard_error(double r, int n) {
    return sqrt((1.0 - r * r) / (n - 2));
}

// Divides a shuffled dataset into smaller groups and calculates the mean
//      for each group. Returns the average of these means and their standard
//      deviation, which gives an idea of how much the results vary between groups
pair<double, double> subsample_mean(const vector<double>& data,
                                    int n_subsamples, mt19937& gen) {
    vector<double> mixed = shuffled_copy(data, gen);
    int subsample_size = mixed.size() / n_subsamples;
    vector<double> means;

    for (int i = 0; i < n_subsamples; ++i) {
        double m = accumulate(mixed.begin() + i * subsample_size,
                              mixed.begin() + (i + 1) * subsample_size,
                              0.0) / subsample_size;
        means.push_back(m);
    }

    return {mean(means), standard_deviation(means)};
}

// Calculates skewness and excess kurtosis for each subsample - 
//      Skewness describes how asymmetric a distribution is, while kurtosis
//      gives an idea of how heavy its tails are compared with a normal distribution

// The function returns the average values and their variation across groups
pair<pair<double, double>, pair<double, double>>
subsample_skew_kurtosis(const vector<double>& data, int n_subsamples,
                        mt19937& gen) {
    vector<double> mixed = shuffled_copy(data, gen);
    int subsample_size = mixed.size() / n_subsamples;
    vector<double> skew_vals, kurt_vals;

    for (int i = 0; i < n_subsamples; ++i) {
        vector<double> sub(mixed.begin() + i * subsample_size,
                           mixed.begin() + (i + 1) * subsample_size);
        double mean_val = mean(sub);

        // Calculate the second, third and fourth central moments
        double m2 = 0.0, m3 = 0.0, m4 = 0.0;

        for (double x : sub) {
            double d = x - mean_val;
            m2 += d * d;
            m3 += d * d * d;
            m4 += d * d * d * d;
        }

        m2 /= sub.size();
        m3 /= sub.size();
        m4 /= sub.size();

        if (m2 == 0.0) continue; // Skewness and kurtosis are undefined here

        skew_vals.push_back(m3 / pow(m2, 1.5));
        kurt_vals.push_back(m4 / (m2 * m2) - 3.0);
    }

    return {
        {mean(skew_vals), standard_deviation(skew_vals)},
        {mean(kurt_vals), standard_deviation(kurt_vals)}
    };
}

// Calculates Pearson's correlation coefficient to measure how closely
//      two sets of measurements follow a linear relationship
double pearson(const vector<double>& x, const vector<double>& y) {
    double mx = mean(x), my = mean(y);
    double sxy = 0.0, sxx = 0.0, syy = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        sxy += (x[i] - mx) * (y[i] - my);
        sxx += (x[i] - mx) * (x[i] - mx);
        syy += (y[i] - my) * (y[i] - my);
    }

    return sxy / sqrt(sxx * syy);
}

// Estimates the correlation and its variation by randomly splitting the
//      paired measurements into groups - the same indices are used for x and y
//      so that each accession stays correctly matched in both datasets
pair<double, double> subsample_pearson(const vector<double>& x,
                                       const vector<double>& y,
                                       int n_subsamples, mt19937& gen) {
    vector<int> idx(x.size());
    iota(idx.begin(), idx.end(), 0);
    shuffle(idx.begin(), idx.end(), gen);

    int subsample_size = idx.size() / n_subsamples;
    vector<double> r_vals;

    for (int i = 0; i < n_subsamples; ++i) {
        vector<double> sx, sy;

        for (int j = i * subsample_size;
             j < (i + 1) * subsample_size; ++j) {
            sx.push_back(x[idx[j]]);
            sy.push_back(y[idx[j]]);
        }

        r_vals.push_back(pearson(sx, sy));
    }

    return {mean(r_vals), standard_deviation(r_vals)};
}


// MAIN FUNCTION - Runs the analysis from start to finish
int main() {
    const int n_subsamples = 10;

    // Using a fixed seed means the random parts of the analysis give
    //      the same results each time the programme is run
    mt19937 gen(42);


    // STEP 1 - Load the measurements for the control and each hormone
    vector<Condition> conditions = {
        {"Control", 359, {}, ""},
        {"IAA",     368, {}, ""},
        {"ABA",     377, {}, ""},
        {"CK",      386, {}, ""}
    };

    cout << "=== STEP 1: Loading data ===" << endl;

    for (auto& c : conditions) {
        string filename = "phenotype_" + to_string(c.phenotype_id) + ".csv";

        if (!download_if_missing(c.phenotype_id, filename)) {
            cerr << "Could not download " << filename
                 << ". Download it manually from "
                 << "https://arapheno.1001genomes.org/rest/phenotype/"
                 << c.phenotype_id << "/values.csv" << endl;
            return 1;
        }

        if (!load_phenotype(filename, c.values, c.pheno_name)) return 1;

        cout << c.label << ": phenotype " << c.phenotype_id
             << " ('" << c.pheno_name << "'), "
             << c.values.size() << " accessions" << endl;
    }

    // Check the printed phenotype names to make sure they refer to
    //        the mean total root length (TRL) measurements
    cout << "(Check that the printed names are the 'Mean(TRL) ...' "
         << "phenotypes, i.e. total root length.)" << endl;


    // STEP 2 - Compare the root length distributions between conditions
    //      Calculates the mean, its uncertainty, skewness and kurtosis
    ofstream summary("results_summary.csv");
    summary << "quantity,condition,N,mean,mean_unc,skewness,skewness_unc,"
            << "kurtosis,kurtosis_unc,sigma_from_1,percent_change,"
            << "mean_full,mean_se,sigma_from_1_se\n";

    cout << fixed << setprecision(3);
    cout << "\n=== STEP 2: Total root length (mm), all accessions per condition ==="
         << endl;

    for (auto& c : conditions) {
        vector<double> v;
        for (auto& kv : c.values) v.push_back(kv.second);

        auto [m, m_unc] = subsample_mean(v, n_subsamples, gen);
        auto [skew_pair, kurt_pair] =
            subsample_skew_kurtosis(v, n_subsamples, gen);

        double m_full = mean(v);         // Mean using every available measurement
        double m_se = standard_error(v); // Standard error using the full dataset

        cout << c.label << " (N = " << v.size() << "): mean = "
             << m << " +/- " << m_unc
             << " mm (SE = " << m_se << "), skewness = "
             << skew_pair.first << " +/- " << skew_pair.second
             << ", kurtosis = " << kurt_pair.first << " +/- "
             << kurt_pair.second << endl;

        summary << "TRL_mm," << c.label << "," << v.size() << ","
                << m << "," << m_unc << ","
                << skew_pair.first << "," << skew_pair.second << ","
                << kurt_pair.first << "," << kurt_pair.second
                << ",,," << m_full << "," << m_se << ",\n";
    }


    // STEP 3 - Compare each accession with itself under the control and
    //          hormone treatments - gives a ratio that accounts for differences
    //          in the starting root lengths of individual accessions

    // A ratio above 1 means longer roots = a ratio below 1 means shorter roots
    const map<int, double>& control = conditions[0].values;
    vector<map<int, double>> ratios(conditions.size());

    // Only calculate a ratio when the same accession has measurements
    //        under both conditions and the control value is greater than zero
    for (size_t h = 1; h < conditions.size(); ++h) {
        for (auto& kv : conditions[h].values) {
            auto it = control.find(kv.first);

            if (it != control.end() && it->second > 0.0) {
                ratios[h][kv.first] = kv.second / it->second;
            }
        }
    }

    cout << "\n=== STEP 3: Relative root length = TRL(hormone) / TRL(control), "
         << "per accession ===" << endl;

    vector<double> ratio_mean(conditions.size(), 0.0),
                   ratio_unc(conditions.size(), 0.0);
    vector<double> ratio_mean_full(conditions.size(), 0.0),
                   ratio_se(conditions.size(), 0.0);

    for (size_t h = 1; h < conditions.size(); ++h) {
        vector<double> r;
        for (auto& kv : ratios[h]) r.push_back(kv.second);

        auto [m, m_unc] = subsample_mean(r, n_subsamples, gen);
        auto [skew_pair, kurt_pair] =
            subsample_skew_kurtosis(r, n_subsamples, gen);

        ratio_mean[h] = m;
        ratio_unc[h] = m_unc;

        // Measure how far the subsampled mean is from 1, the value
        //      expected if the hormone has no effect on root length
        double sigma_from_1 = (m - 1.0) / m_unc;
        double percent_change = (m - 1.0) * 100.0;

        // Also calculate the mean and standard error using every ratio
        double r_mean_full = mean(r);
        double r_se = standard_error(r);
        double sigma_se = (r_mean_full - 1.0) / r_se;

        ratio_mean_full[h] = r_mean_full;
        ratio_se[h] = r_se;

        cout << conditions[h].label << " (N = " << r.size()
             << "): mean ratio = " << m << " +/- " << m_unc
             << " (SE = " << r_se << ", "
             << showpos << sigma_se << noshowpos << " SE from 1)"
             << "  (" << showpos << percent_change << "%"
             << noshowpos << ")"
             << ", deviation from 1 = "
             << showpos << sigma_from_1 << " sigma" << noshowpos
             << (fabs(sigma_from_1) > 3
                 ? "  -> significant effect (>3 sigma)"
                 : "  -> consistent with no effect")
             << endl;

        cout << "    skewness = " << skew_pair.first << " +/- "
             << skew_pair.second << ", kurtosis = "
             << kurt_pair.first << " +/- " << kurt_pair.second << endl;

        summary << "ratio_to_control," << conditions[h].label << ","
                << r.size() << "," << m << "," << m_unc << ","
                << skew_pair.first << "," << skew_pair.second << ","
                << kurt_pair.first << "," << kurt_pair.second << ","
                << sigma_from_1 << "," << percent_change << ","
                << r_mean_full << "," << r_se << "," << sigma_se << "\n";
    }

    summary.close();


    // STEP 4 - Compare the average ratios for each pair of hormones

    // The combined uncertainty is calculated from the two individual
    //      uncertainties - a difference greater than 3 times this uncertainty
    //      is flagged by the programme as a potentially meaningful difference
    ofstream pairwise("pairwise_comparisons.csv");
    pairwise << "hormone_1,hormone_2,mean_difference,combined_unc,"
             << "sigma,diff_full,combined_se,sigma_se\n";

    cout << "\n=== STEP 4: Difference between hormones (3 sigma criterion) ==="
         << endl;

    for (size_t a = 1; a < conditions.size(); ++a) {
        for (size_t b = a + 1; b < conditions.size(); ++b) {
            double diff = ratio_mean[a] - ratio_mean[b];

            // Combine the uncertainties, assuming they are independent
            double combined = sqrt(ratio_unc[a] * ratio_unc[a] +
                                   ratio_unc[b] * ratio_unc[b]);
            double sigma = fabs(diff) / combined;

            // Repeat the comparison using the means and standard errors
            //      calculated from all available ratios
            double diff_full = ratio_mean_full[a] - ratio_mean_full[b];
            double combined_se = sqrt(ratio_se[a] * ratio_se[a] +
                                      ratio_se[b] * ratio_se[b]);
            double sigma_se = fabs(diff_full) / combined_se;

            cout << conditions[a].label << " vs "
                 << conditions[b].label << ": difference = "
                 << diff << " +/- " << combined << " ("
                 << sigma << " sigma) -> "
                 << (sigma > 3
                     ? "statistically different"
                     : "statistically consistent")
                 << endl;

            pairwise << conditions[a].label << ","
                     << conditions[b].label << ","
                     << diff << "," << combined << "," << sigma << ","
                     << diff_full << "," << combined_se << ","
                     << sigma_se << "\n";
        }
    }

    pairwise.close();


    // STEP 5 - Check whether accessions that respond strongly to one
    //          hormone also tend to respond strongly to another

    // Pearson's correlation is calculated using the same accessions
    //          for both hormones, so the pairs remain matched
    ofstream corr_file("correlations.csv");
    corr_file << "hormone_1,hormone_2,N,pearson_r,pearson_r_unc,"
              << "pearson_r_full,pearson_r_se\n";

    cout << "\n=== STEP 5: Correlation of the responses between hormones "
         << "(same accessions) ===" << endl;

    for (size_t a = 1; a < conditions.size(); ++a) {
        for (size_t b = a + 1; b < conditions.size(); ++b) {
            vector<double> x, y;

            for (auto& kv : ratios[a]) {
                auto it = ratios[b].find(kv.first);

                // Keep only accessions with a ratio for both hormones
                if (it != ratios[b].end()) {
                    x.push_back(kv.second);
                    y.push_back(it->second);
                }
            }

            // There must be enough paired measurements to from the groups
            if (x.size() < (size_t)(2 * n_subsamples)) continue;

            auto [r, r_unc] = subsample_pearson(x, y, n_subsamples, gen);

            // Calculate the correlation and its standard error using
            //          all the paired measurements as well
            double r_full = pearson(x, y);
            double r_se = pearson_standard_error(r_full, x.size());

            cout << conditions[a].label << " vs "
                 << conditions[b].label << " (N = " << x.size()
                 << "): r = " << r << " +/- " << r_unc
                 << " (all pairs: r = " << r_full
                 << ", SE = " << r_se << ")" << endl;

            corr_file << conditions[a].label << ","
                      << conditions[b].label << "," << x.size() << ","
                      << r << "," << r_unc << ","
                      << r_full << "," << r_se << "\n";
        }
    }

    corr_file.close();


    // STEP 6 - Save the measurements for individual accessions and prepare
    //          a histogram table showing how the hormone-to-control ratios are spread - 
    //          These files can be opened in a spreadsheet programme to make graphs! :)
    ofstream per_acc("ratios_per_accession.csv");
    per_acc << "accession_id,TRL_control,TRL_IAA,TRL_ABA,TRL_CK,"
            << "ratio_IAA,ratio_ABA,ratio_CK\n";

    for (auto& kv : control) {
        int id = kv.first;
        per_acc << id << "," << kv.second;

        // Save the original root lengths. Leave the cell empty if a
        // measurement is missing for that accession and treatment.
        for (size_t h = 1; h < conditions.size(); ++h) {
            auto it = conditions[h].values.find(id);
            per_acc << ","
                    << (it != conditions[h].values.end()
                        ? to_string(it->second) : "");
        }

        // Save the hormone-to-control ratios in the same way
        for (size_t h = 1; h < conditions.size(); ++h) {
            auto it = ratios[h].find(id);
            per_acc << ","
                    << (it != ratios[h].end()
                        ? to_string(it->second) : "");
        }

        per_acc << "\n";
    }

    per_acc.close();

    // Group the ratios into 20 bins between 0 and 2
    //          Values below 0 or above 2 are included in the first or last bin
    const int n_bins = 20;
    const double lo = 0.0, hi = 2.0;
    const double bin_width = (hi - lo) / n_bins;

    vector<vector<int>> hist(conditions.size(), vector<int>(n_bins, 0));

    for (size_t h = 1; h < conditions.size(); ++h) {
        for (auto& kv : ratios[h]) {
            int bin = (int)floor((kv.second - lo) / bin_width);
            bin = max(0, min(n_bins - 1, bin));
            hist[h][bin]++;
        }
    }

    // Write the number of accessions in each bin to a CSV file
    ofstream hist_file("ratio_histogram.csv");
    hist_file << "bin_low,bin_high,IAA,ABA,CK\n";

    for (int b = 0; b < n_bins; ++b) {
        hist_file << lo + b * bin_width << ","
                  << lo + (b + 1) * bin_width << ","
                  << hist[1][b] << ","
                  << hist[2][b] << ","
                  << hist[3][b] << "\n";
    }

    hist_file.close();

    // Print the names of all the files created by the analyssi
    cout << "\nFiles written: results_summary.csv, pairwise_comparisons.csv, "
         << "correlations.csv, ratios_per_accession.csv, "
         << "ratio_histogram.csv" << endl;

    return 0;
}
