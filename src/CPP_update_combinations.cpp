
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void updateDistancesWithCombinations(NumericVector& length_root_distances,
                                                                             NumericVector& topological_root_distances,
                                                                             IntegerVector& left_partition,
                                                                             IntegerVector& right_partition,
                                                                             IntegerVector& index_offsets,
                                                                             double distance_to_root,
                                                                             int edges_to_root)
{
  // Iterate through all combinations.
  for(int i=0; i < left_partition.size(); ++i) {
    for(int j=0; j < right_partition.size(); ++j) {
      int first_leaf = left_partition[i];
      int second_leaf = right_partition[j];
      // Because of the symmetric distances.
      if(left_partition[i] > right_partition[j]) {
        first_leaf = right_partition[j];
        second_leaf = left_partition[i];
      }
      // Roll the index (notice we take into account C++ indices here, starting at 0).
      int combination_index = index_offsets[first_leaf-1] + (second_leaf - first_leaf) - 1;
      // Update the vectors.
      length_root_distances[combination_index] = distance_to_root;
      topological_root_distances[combination_index] = edges_to_root;
    }
  }
}