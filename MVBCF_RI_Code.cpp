//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppDist.h>


using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]

// Node class definition
class Node {
  
public:
  
  //Attributes
  arma::colvec mu;
  int variable;
  double split_val;
  arma::uvec observations;
  arma::uvec test_observations;
  bool is_terminal;
  bool in_use;
  
  // Default Constructor
  Node() {
    variable = -1;
    split_val = -1;
    is_terminal = false;
    in_use = false;
  }
  
  // Copy constructor for the Node class
  Node(const Node& other) {
    mu = other.mu;
    variable = other.variable;
    split_val = other.split_val;
    observations = other.observations;
    test_observations = other.test_observations;
    is_terminal = other.is_terminal;
    in_use = other.in_use;
  }
  
  
  // Method for updating mu
  void update_mu(arma::mat sigma, 
                 arma::mat sigma_mu, 
                 arma::mat y_resid) {
    
    
    arma::colvec mu_0(y_resid.n_cols, arma::fill::zeros);
    
    double nj = sum(observations);
    
    arma::mat node_resid = y_resid.rows(find(observations == 1));
    
    arma::colvec y_bar = arma::mean(node_resid, 0).t();
    
    arma::mat part1 = arma::inv(arma::inv(sigma_mu) + nj * arma::inv(sigma));
    
    arma::vec part2 = arma::inv(sigma_mu) * mu_0 + nj * arma::inv(sigma) * y_bar;
    
    arma::rowvec temp = rmvnorm(1, part1 * part2, part1);
      
    mu = temp.t();
      
  }
  
  
  // Method for updating tau
  void update_tau(arma::mat sigma, 
                 arma::mat sigma_tau, 
                 arma::mat y_resid,
                 arma::mat Z) {
    
    
    arma::colvec mu_0(y_resid.n_cols, arma::fill::zeros);
    
    arma::mat node_z = Z.rows(find(observations == 1));
    
    arma::mat njz = node_z.t()*node_z;
    
    arma::mat node_resid = y_resid.rows(find(observations == 1));
    
    arma::mat part1 = arma::inv( arma::inv(sigma_tau) + njz % arma::inv(sigma) );
    
    arma::mat tricky_bit(y_resid.n_cols, 1, arma::fill::zeros);
    
    arma::mat node_z_t = node_z.t();
    
    arma::mat node_resid_t = node_resid.t();
    
    for(int i = 0; i<node_resid.n_rows; i++)
    {
      tricky_bit = tricky_bit + node_z_t.col(i) % (arma::inv(sigma)*node_resid_t.col(i));
    }
    
    arma::mat part2 = arma::inv(sigma_tau) * mu_0 + tricky_bit;
    
    arma::rowvec temp = rmvnorm(1, part1 * part2, part1);
    
    mu = temp.t();
    
  }
  
  
};





// Tree class definition
class Tree {
  
public:
  
  std::vector<Node> node_vector;
  
  // Constructor
  Tree(int num_nodes = 1, int num_obs = 1, int num_test_obs = 1) {
    node_vector.resize(num_nodes);
    node_vector[0].observations = arma::uvec(num_obs, arma::fill::ones);
    node_vector[0].test_observations = arma::uvec(num_test_obs, arma::fill::ones);
    node_vector[0].in_use=true;
    node_vector[0].is_terminal=true;
  }
  
  // Copy constructor for the Tree class
  Tree(const Tree& other) {
    // Copy each node from the other tree
    for (const Node& node : other.node_vector) {
      node_vector.push_back(Node(node));  // Invoke the Node copy constructor
    }
  }

  
  // Method for updating all terminal nodes
  void update_nodes(arma::mat sigma, 
                    arma::mat sigma_mu, 
                    arma::mat y_resid) {
    
    int num_nodes = node_vector.size();
    
    for(int i=0; i<num_nodes; i++)
    {
      if(node_vector[i].is_terminal & node_vector[i].in_use)
      {
        node_vector[i].update_mu(sigma, sigma_mu, y_resid);
      }
    }
  }
  
  
  // Method for updating all terminal nodes
  void update_nodes_tau(arma::mat sigma, 
                    arma::mat sigma_tau, 
                    arma::mat y_resid,
                    arma::mat Z) {
    
    int num_nodes = node_vector.size();
    
    for(int i=0; i<num_nodes; i++)
    {
      if(node_vector[i].is_terminal & node_vector[i].in_use)
      {
        node_vector[i].update_tau(sigma, sigma_tau, y_resid, Z);
      }
    }
  }

  
  // Method for selecting a terminal node
  int get_terminal_node() {
    
    std::vector<int> valid_indices;
    
    for (int i = 0; i < node_vector.size(); i++) {
      if (node_vector[i].is_terminal & node_vector[i].in_use) {
        valid_indices.push_back(i);
      }
    }
    
    return valid_indices[floor(R::runif(0, valid_indices.size()))];
    
  }
  
  
  // Method for selecting a non terminal node
  int get_non_terminal_node() {
    
    std::vector<int> valid_indices;
    
    for (int i = 0; i < node_vector.size(); i++) {
      if (!node_vector[i].is_terminal & node_vector[i].in_use) {
        valid_indices.push_back(i);
      }
    }
    
    if(valid_indices.size()>0)
    {
      return valid_indices[floor(R::runif(0, valid_indices.size()))];
    }
    else
    {
      return -1;
    }
    
  }
  
  
  // Method for selecting a non terminal node with a parent
  int get_parent_child() {
    
    std::vector<int> valid_indices;
    
    for (int i = 0; i < node_vector.size(); i++) {
      if (!node_vector[i].is_terminal & node_vector[i].in_use & i!=0) {
        valid_indices.push_back(i);
      }
    }
    
    if(valid_indices.size()>0)
    {
      return valid_indices[floor(R::runif(0, valid_indices.size()))];
    }
    else
    {
      return -1;
    }
    
  }
  
  
  // Method for selecting a parent of two terminal nodes
  int get_terminal_parent() {
    
    std::vector<int> valid_indices;
    
    for (int i = 0; i < node_vector.size(); i++) {
      if(node_vector[i].in_use & !node_vector[i].is_terminal)
      {
        if(node_vector[2*i+1].in_use & node_vector[2*i+1].is_terminal & node_vector[2*i+2].in_use & node_vector[2*i+2].is_terminal)
        {
          valid_indices.push_back(i);
        }
      }
    }
    
    if(valid_indices.size()>0)
    {
      return valid_indices[floor(R::runif(0, valid_indices.size()))];
    }
    else
    {
      return -1;
    }
  }
  
  
  
  // Method for growing tree
  void grow(arma::mat X, arma::mat X_test, int p, int min_nodesize) {
    
    int grow_index = get_terminal_node();
    
    int variable = floor(R::runif(0, p));
    
    node_vector[grow_index].variable = variable;
    
    arma::colvec X_col = X.col(variable);
    arma::colvec X_test_col = X_test.col(variable);
    
    arma::colvec X_col_subset = X_col.rows(find(node_vector[grow_index].observations == 1));
    arma::colvec X_test_col_subset = X_test_col.rows(find(node_vector[grow_index].test_observations == 1));
    
    arma::colvec X_unique = arma::unique(X_col_subset);
    
    double split_val;
    
    if(X_unique.n_rows>0)
    {
      int random_index = floor(R::runif(0, X_unique.n_rows));
      
      split_val = X_unique(random_index);
    }
    else
    {
      split_val = -1;
    }
    
    node_vector[grow_index].split_val = split_val;
    
    arma::uvec is_less = X_col<=split_val;
    arma::uvec is_less_test = X_test_col<=split_val;
    
    arma::uvec less_subset = node_vector[grow_index].observations && is_less;
    
    arma::uvec more_subset = node_vector[grow_index].observations && (1-is_less);
    
    int sum_less = sum(less_subset);
    
    int sum_more = sum(more_subset);
    
    if(sum_more>=min_nodesize & sum_less>=min_nodesize)
    {
      if(node_vector.size()<2*grow_index+2+1)
      {
        node_vector.resize(2*grow_index+2+1);
      }
      
      int child_left = 2*grow_index+1;
      
      int child_right = 2*grow_index+2;
      
      node_vector[child_left].observations = node_vector[grow_index].observations && is_less;
      node_vector[child_left].test_observations = node_vector[grow_index].test_observations && is_less_test;
      node_vector[child_left].is_terminal = true;
      node_vector[child_left].in_use = true;
      
      node_vector[child_right].observations = node_vector[grow_index].observations && (1-is_less);
      node_vector[child_right].test_observations = node_vector[grow_index].test_observations && (1-is_less_test);
      node_vector[child_right].is_terminal = true;
      node_vector[child_right].in_use = true;
      
      node_vector[grow_index].is_terminal = false;
      node_vector[grow_index].in_use = true;
    }
  }
  
  
  // Method for pruning
  void prune() {
    
    int prune_index = get_terminal_parent();
    
    if(prune_index!=-1)
    {
      node_vector[prune_index*2+1].in_use = false;
      node_vector[prune_index*2+1].is_terminal = false;
      
      node_vector[prune_index*2+2].in_use = false;
      node_vector[prune_index*2+2].is_terminal = false;
      
      node_vector[prune_index].is_terminal = true;
      node_vector[prune_index].in_use = true;
    }
  }
  
  
  
  // Method for updating observations
  void change_update(arma::mat X, arma::mat X_test) {
    
    int num_nodes = node_vector.size();
    
    for(int i = 0; i < num_nodes; i++)
    {
      if(!node_vector[i].is_terminal & node_vector[i].in_use)
      {
        int child_left = 2*i+1;
        int child_right = 2*i+2;
        
        int variable = node_vector[i].variable;
        double split_val = node_vector[i].split_val;
        
        arma::uvec is_less = X.col(variable)<=split_val;
        arma::uvec is_more = X.col(variable)>split_val;
        
        arma::uvec is_less_test = X_test.col(variable)<=split_val;
        arma::uvec is_more_test = X_test.col(variable)>split_val;
        
        node_vector[child_left].observations = node_vector[i].observations && is_less;
        node_vector[child_right].observations = node_vector[i].observations && is_more;
        
        node_vector[child_left].test_observations = node_vector[i].test_observations && is_less_test;
        node_vector[child_right].test_observations = node_vector[i].test_observations && is_more_test;
      }
    }
  }
  
  
  // Method for changing
  void change(arma::mat X, int p) {
    
    int change_index = get_non_terminal_node();
    
    if(change_index!=-1)
    {
      int variable = floor(R::runif(0, p));
      
      node_vector[change_index].variable = variable;
      
      arma::colvec X_col = X.col(variable);
      
      X_col = X_col.rows(find(node_vector[change_index].observations == 1));
      
      arma::colvec X_unique = arma::unique(X_col);
      
      if(X_unique.size()>0)
      {
        int random_index = floor(R::runif(0, X_unique.n_rows));
        
        node_vector[change_index].split_val = X_unique(random_index);
      }
      else
      {
        node_vector[change_index].split_val = -1;
      }
    }
  }

  
  
  // Method for swapping
  void swap() {
    
    int swap_index = get_parent_child();
    
    if(swap_index!=-1)
    {
      int parent_index = (swap_index-1)/2;
      
      int parent_variable = node_vector[parent_index].variable;
      double parent_split_val = node_vector[parent_index].split_val;
      
      int child_variable = node_vector[swap_index].variable;
      double child_split_val = node_vector[swap_index].split_val;
      
      node_vector[parent_index].variable = child_variable;
      node_vector[parent_index].split_val = child_split_val;
      
      node_vector[swap_index].variable = parent_variable;
      node_vector[swap_index].split_val = parent_split_val;
    }
  }
  
  
  // Method for checking if any nodes are empty
  bool has_empty_nodes(int min_nodesize) {
    
    int num_nodes = node_vector.size();
    
    for(int i=0; i<num_nodes; i++)
    {
      if(node_vector[i].in_use & node_vector[i].is_terminal)
      {
        if(sum(node_vector[i].observations)<min_nodesize)
        {
          return true;
        } 
      }
    }
    
    return false;
  }

  double log_lik(arma::mat sigma_mu, 
                 arma::mat sigma,
                 double alpha,
                 double beta,
                 arma::mat y_resid){
    
    
    double log_lik = 0.0;
    
    for(int i = 0; i < node_vector.size(); i++)
    {
      if(node_vector[i].in_use & node_vector[i].is_terminal)
      {
        double nj = sum(node_vector[i].observations);
        arma::mat node_resid = y_resid.rows(find(node_vector[i].observations == 1));
        
        arma::mat sigma_j0_inv = nj*arma::inv(sigma) + arma::inv(sigma_mu);
        arma::mat mu_j0 = (arma::inv(sigma_j0_inv))*(arma::inv(sigma))*(arma::sum(node_resid, 0).as_col());
        
        double eq1p1 = (-1.0*nj/2.0)*log(arma::det(sigma));
        double eq1p2 = (-1.0/2.0)*log(arma::det(sigma_mu));
        double eq1p3 = (-1.0/2.0)*log(arma::det(arma::inv(sigma_mu)+nj*arma::inv(sigma)));
        double eq1p4 = (-1.0/2.0)*arma::accu((node_resid.t()*node_resid)%arma::inv(sigma));
        
        double eq1p5 = arma::accu((1.0/2.0)*(mu_j0.t())*(sigma_j0_inv)*(mu_j0));
        
        double eq4p1 = log(1.0-alpha*pow(1+floor(log2(i + 1)), (-1*beta)));
        
        log_lik += eq1p1 + eq1p2 + eq1p3 + eq1p4 + eq1p5 + eq4p1;
      }
      else if(node_vector[i].in_use & !node_vector[i].is_terminal)
      {
        double eq4p2 = log(alpha)-beta*log(1+floor(log2(i + 1)));
        
        log_lik += eq4p2;
      }
    }
    
    return log_lik;
  }
  
  
  double log_lik_tau(arma::mat sigma_tau, 
                 arma::mat sigma,
                 double alpha,
                 double beta,
                 arma::mat y_resid,
                 arma::mat Z){
    
    
    double log_lik = 0.0;
    
    for(int i = 0; i < node_vector.size(); i++)
    {
      if(node_vector[i].in_use & node_vector[i].is_terminal)
      {
        double nj = sum(node_vector[i].observations);
        arma::mat node_resid = y_resid.rows(find(node_vector[i].observations == 1));
        arma::mat node_z = Z.rows(find(node_vector[i].observations == 1));
        
        arma::mat tricky_bit(y_resid.n_cols, 1, arma::fill::zeros);
        
        arma::mat node_z_t = node_z.t();
        
        arma::mat node_resid_t = node_resid.t();
        
        for(int t = 0; t<node_resid.n_rows; t++)
        {
          tricky_bit = tricky_bit + node_z_t.col(t) % (arma::inv(sigma)*node_resid_t.col(t));
        }
        
        arma::mat sigma_j0_inv = (node_z.t()*node_z)%arma::inv(sigma) + arma::inv(sigma_tau);
        arma::mat tau_j0 = (arma::inv(sigma_j0_inv))*tricky_bit;
        
        double eq1p1 = (-1.0*nj/2.0)*log(arma::det(sigma));
        double eq1p2 = (-1.0/2.0)*log(arma::det(sigma_tau));
        double eq1p3 = (-1.0/2.0)*log(arma::det(arma::inv(sigma_tau)+(node_z.t() * node_z)%arma::inv(sigma)));
        double eq1p4 = (-1.0/2.0)*arma::accu((node_resid.t()*node_resid)%arma::inv(sigma));
        
        double eq1p5 = arma::accu((1.0/2.0)*(tau_j0.t())*(sigma_j0_inv)*(tau_j0));
        
        double eq4p1 = log(1.0-alpha*pow(1+floor(log2(i + 1)), (-1*beta)));
        
        log_lik += eq1p1 + eq1p2 + eq1p3 + eq1p4 + eq1p5 + eq4p1;
      }
      else if(node_vector[i].in_use & !node_vector[i].is_terminal)
      {
        double eq4p2 = log(alpha)-beta*log(1+floor(log2(i + 1)));
        
        log_lik += eq4p2;
      }
    }
    
    return log_lik;
  }
  
  


  
  // Method for getting predictions from tree
  arma::mat get_predictions(int num_outcomes) {
    int num_obs = node_vector[0].observations.n_elem;
    int num_nodes = node_vector.size();
    
    arma::mat predictions(num_obs, num_outcomes, arma::fill::zeros);
    
    for(int i = 0; i < num_nodes; i++) {
      for(int j = 0; j < num_obs; j++) {
        if(node_vector[i].is_terminal & node_vector[i].in_use && (node_vector[i].observations[j]==1)) {
          for(int k=0; k<num_outcomes; k++)
          {
            predictions(j, k) = node_vector[i].mu(k);
          }
        }
      }
    }
    
    return predictions;
  }
  
  arma::mat get_test_predictions(int num_outcomes) {
    int num_obs = node_vector[0].test_observations.n_elem;
    int num_nodes = node_vector.size();
    
    arma::mat test_predictions(num_obs, num_outcomes, arma::fill::zeros);
    
    for(int i = 0; i < num_nodes; i++) {
      for(int j = 0; j < num_obs; j++) {
        if(node_vector[i].is_terminal & node_vector[i].in_use && (node_vector[i].test_observations[j]==1)) {
          for(int k=0; k<num_outcomes; k++)
          {
            test_predictions(j, k) = node_vector[i].mu(k);
          }
        }
      }
    }
    
    return test_predictions;
  }
  
};


// Forest class definition
class Forest {
  
public:
  
  std::vector<Tree> tree_vector;
  
  // Constructor
  Forest(int num_trees=1, int num_nodes = 1, int num_obs=1, int num_test_obs = 1) {
    
    tree_vector.resize(num_trees);
    
    for(int i=0; i<num_trees; i++)
    {
      tree_vector[i] = Tree(num_nodes, num_obs, num_test_obs);
    }
  }
};


arma::mat sum_over_cube_without_slice(arma::cube arma_cube, int slice_removed) {
  // Get the dimensions of the cube
  int n_rows = arma_cube.n_rows;
  int n_cols = arma_cube.n_cols;
  int n_layers = arma_cube.n_slices;
  
  // Create a matrix to store the result
  arma::mat result_matrix(n_rows, n_cols, arma::fill::zeros);
  
  // Sum over the third dimension (layers)
  for (int i = 0; i < n_layers; i++) {
    if(i != slice_removed)
    {
      result_matrix += arma_cube.slice(i); 
    }
  }
 
 return result_matrix;
}

 

arma::mat sample_sigma(double n,
                    int v_0,
                    arma::mat y,
                    arma::mat preds,
                    arma::mat sigma_0) {
  
  arma::mat resid = y - preds;
  arma::mat sig_mat = resid.t()*resid;
  
  return riwish(v_0+n, sigma_0+sig_mat);
}


//Functions for the hierarchical parameters:
arma::colvec sample_alpha(arma::mat t,
                          arma::colvec m,
                          arma::mat s,
                          arma::mat alpha_is)
{
  int G = alpha_is.n_rows;
  
  arma::colvec alpha_bar = mean(alpha_is, 0).t();
  
  arma::mat part1 = arma::inv(arma::inv(t) + G*arma::inv(s));
  
  arma::mat part2 = arma::inv(t)*m + G*arma::inv(s)*alpha_bar;
  
  arma::rowvec temp = rmvnorm(1, part1 * part2, part1);
  
  arma::colvec alpha = temp.t();
  
  return alpha;
  
}

arma::mat sample_alpha_is(arma::mat t,
                          arma::mat y_vals,
                          arma::mat sigma,
                          arma::colvec alpha_mean,
                          arma::rowvec group_ids)
{
  arma::rowvec unique_elements = arma::unique(group_ids);
  
  int G = unique_elements.n_elem;
  
  arma::mat alpha_is(G, y_vals.n_cols);
  
  for(int i = 0; i<G; i++)
  {
    arma::mat y_vals_g = y_vals.rows(find(group_ids == i));
    
    int ng = y_vals_g.n_rows;
    
    arma::colvec y_bar = mean(y_vals_g, 0).t();
    
    arma::mat part1 = arma::inv(arma::inv(sigma) + ng*arma::inv(t));
    
    arma::mat part2 = arma::inv(sigma)*alpha_mean + ng*arma::inv(t)*y_bar;
    
    arma::rowvec alpha_g = rmvnorm(1, part1*part2, part1);
    
    alpha_is.row(i) = alpha_g;
  }
  
  return alpha_is;
}

arma::mat sample_t(double b,
                   arma::mat alpha_is,
                   arma::colvec alpha_mean,
                   arma::mat Sb)
{
  int G = alpha_is.n_rows;
  
  arma::mat alpha_resid = alpha_is.each_row() - alpha_mean.t();
  
  arma::mat Sg = alpha_resid.t()*alpha_resid;
  
  arma::mat temp = riwish(b+G, Sb + Sg);
  
  return temp;
}


// [[Rcpp::export]]
List fast_bart(arma::mat X_con,
               arma::mat y,
               arma::mat Z,
               arma::mat X_mod,
               arma::mat X_con_test,
               arma::mat X_mod_test,
               double alpha,
               double beta,
               double alpha_tau,
               double beta_tau,
               arma::mat sigma_mu,
               arma::mat sigma_tau,
               int v_0,
               arma::mat sigma_0,
               int n_iter,
               int n_tree,
               int n_tree_tau,
               int min_nodesize,
               arma::rowvec group_ids,
               arma::rowvec group_ids_test,
               arma::mat Sb,
               arma::colvec m,
               arma::mat s,
               double b)
{
  
  
  
  //set sigma
  arma::mat sigma(y.n_cols, y.n_cols, arma::fill::eye);
  
  
  //scale y
  arma::rowvec col_means = mean(y, 0);
  arma::rowvec col_stdev = stddev(y, 0);
  arma::mat y_scaled = y.each_row() - col_means;
  y_scaled.each_row() /= col_stdev;
  
  //get number of variables p, and rows n, and outcomes 
  int n = y_scaled.n_rows;
  int n_test = X_con_test.n_rows;
  int p = X_con.n_cols;
  int p_tau = X_mod.n_cols;
  int outcomes = y_scaled.n_cols;
   
  
  
  
  //For holding tree predictions at each iteration
  arma::cube tree_preds(n, y.n_cols, n_tree);
  arma::cube tree_preds_tau(n, y.n_cols, n_tree_tau);
  
  arma::cube tree_preds_test(n_test, y.n_cols, n_tree);
  arma::cube tree_preds_tau_test(n_test, y.n_cols, n_tree_tau);
  
  //For holding overall predictions from each iteration
  arma::cube preds_mat(n, y.n_cols, n_iter);
  arma::cube preds_mat_tau(n, y.n_cols, n_iter);
  
  arma::cube preds_mat_test(n_test, y.n_cols, n_iter);
  arma::cube preds_mat_tau_test(n_test, y.n_cols, n_iter);
  
  //For holding tau from each iteration
  arma::cube sigmas(y.n_cols, y.n_cols, n_iter);
  
  StringVector choices = {"Grow", "Prune", "Change", "Swap"};
  
  
  //initialise the forests
  Forest bart_forest(n_tree, 1, n, n_test);
  
  Forest tau_forest(n_tree_tau, 1, n, n_test);
  
  
  //Stuff related to hierarchical parameters
  arma::rowvec unique_elements = arma::unique(group_ids);
  int num_groups = unique_elements.n_elem;
  
  //initialise updated parameters
  arma::mat t(outcomes, outcomes, arma::fill::eye);
  arma::colvec alpha_mean(outcomes, arma::fill::zeros);
  arma::mat alpha_is(num_groups, outcomes, arma::fill::zeros);
  arma::mat alpha_preds(n, outcomes, arma::fill::zeros);
  arma::mat alpha_preds_test(n_test, outcomes, arma::fill::zeros);
  
  
  //set up storage
  arma::cube ts(outcomes, outcomes, n_iter, arma::fill::zeros);
  arma::mat alpha_means(outcomes, n_iter, arma::fill::zeros);
  arma::cube random_effects(num_groups, outcomes, n_iter, arma::fill::zeros);
  
  
  
  
  for(int iter = 0; iter < n_iter; iter++)
  {
    for(int tree_num = 0; tree_num < n_tree; tree_num++)
    {
      arma::mat y_resid = y_scaled-sum_over_cube_without_slice(tree_preds, tree_num)
                                  -Z%sum_over_cube_without_slice(tree_preds_tau, -1)
                                  -alpha_preds;
      
      String choice = sample(choices, 1)[0];
      
      Tree proposal_tree = Tree(bart_forest.tree_vector[tree_num]);
      
      
      if(choice == "Grow")
      {
        proposal_tree.grow(X_con, X_con_test, p, min_nodesize);
      }
      
      
      if(choice == "Prune")
      {
        proposal_tree.prune();
      }
      
      
      if(choice == "Change")
      {
        proposal_tree.change(X_con, p);
        proposal_tree.change_update(X_con, X_con_test);
      }
      
      
      if(choice == "Swap")
      {
        proposal_tree.swap();
        proposal_tree.change_update(X_con, X_con_test);
      }
      
      
      if(!proposal_tree.has_empty_nodes(min_nodesize))
      {
        double lnew = proposal_tree.log_lik(sigma_mu, 
                                            sigma,
                                            alpha,
                                            beta,
                                            y_resid);
          
        double lold = bart_forest.tree_vector[tree_num].log_lik(sigma_mu, 
                                                                sigma,
                                                                alpha,
                                                                beta,
                                                                y_resid);
          
        double a = exp(lnew-lold);
        if(a > R::runif(0, 1))
        {
          bart_forest.tree_vector[tree_num] = Tree(proposal_tree);
        }
      }
      
      
      
      bart_forest.tree_vector[tree_num].update_nodes(sigma, sigma_mu, y_resid);
      
      arma::mat tree_preds_from_iter = bart_forest.tree_vector[tree_num].get_predictions(outcomes);
      arma::mat tree_preds_from_iter_test = bart_forest.tree_vector[tree_num].get_test_predictions(outcomes);
      
      tree_preds.slice(tree_num) = tree_preds_from_iter;
      tree_preds_test.slice(tree_num) = tree_preds_from_iter_test;
      
    
    }
    for(int tree_num = 0; tree_num < n_tree_tau; tree_num++)
    {
      arma::mat y_resid = y_scaled-sum_over_cube_without_slice(tree_preds, -1)
                                  -Z%sum_over_cube_without_slice(tree_preds_tau, tree_num)
                                  -alpha_preds;
      
      String choice = sample(choices, 1)[0];
      
      Tree proposal_tree = Tree(tau_forest.tree_vector[tree_num]);
      
      
      if(choice == "Grow")
      {
        proposal_tree.grow(X_mod, X_mod_test, p_tau, min_nodesize);
      }
      
      
      if(choice == "Prune")
      {
        proposal_tree.prune();
      }
      
      
      if(choice == "Change")
      {
        proposal_tree.change(X_mod, p_tau);
        proposal_tree.change_update(X_mod, X_mod_test);
      }
      
      
      if(choice == "Swap")
      {
        proposal_tree.swap();
        proposal_tree.change_update(X_mod, X_mod_test);
      }
      
      
      if(!proposal_tree.has_empty_nodes(min_nodesize))
      {
        double lnew = proposal_tree.log_lik_tau(sigma_tau, 
                                            sigma,
                                            alpha_tau,
                                            beta_tau,
                                            y_resid,
                                            Z);
        
        double lold = tau_forest.tree_vector[tree_num].log_lik_tau(sigma_tau, 
                                                                sigma,
                                                                alpha_tau,
                                                                beta_tau,
                                                                y_resid,
                                                                Z);
        
        double a = exp(lnew-lold);
        if(a > R::runif(0, 1))
        {
          tau_forest.tree_vector[tree_num] = Tree(proposal_tree);
        }
      }
      
      
      
      tau_forest.tree_vector[tree_num].update_nodes_tau(sigma, sigma_tau, y_resid, Z);
      
      arma::mat tree_preds_from_iter = tau_forest.tree_vector[tree_num].get_predictions(outcomes);
      arma::mat tree_preds_from_iter_test = tau_forest.tree_vector[tree_num].get_test_predictions(outcomes);
      
      tree_preds_tau.slice(tree_num) = tree_preds_from_iter;
      tree_preds_tau_test.slice(tree_num) = tree_preds_from_iter_test;
      
      
    }
      
    
    Rcpp::Rcout << "Total of " << iter+1 << " of " << n_iter << " iterations completed! " << "(" << (float)(iter+1)/(float)n_iter*100 << "%)             " << "\r";
    Rcpp::Rcout.flush();
    
    arma::mat iter_preds = sum_over_cube_without_slice(tree_preds, -1);
    arma::mat iter_preds_tau = sum_over_cube_without_slice(tree_preds_tau, -1);
    
    arma::mat iter_preds_test = sum_over_cube_without_slice(tree_preds_test, -1);
    arma::mat iter_preds_tau_test = sum_over_cube_without_slice(tree_preds_tau_test, -1);
  
    arma::mat y_resid = y_scaled - iter_preds - Z%iter_preds_tau;
    
    alpha_is = sample_alpha_is(t,
                               y_resid,
                               sigma,
                               alpha_mean,
                               group_ids);
    
    alpha_mean = sample_alpha(t,
                              m,
                              s,
                              alpha_is);
    
    
    t = sample_t(b,
                 alpha_is,
                 alpha_mean,
                 Sb);
    
    for(int i = 0; i<n; i++)
    {
      arma::rowvec alpha_row = alpha_is.row(group_ids(i));
      
      alpha_preds.row(i) = alpha_row;
    }
    for(int i = 0; i<n_test; i++)
    {
      arma::rowvec alpha_row_test = alpha_is.row(group_ids_test(i));
      
      alpha_preds_test.row(i) = alpha_row_test;
    }
    
    iter_preds = iter_preds+alpha_preds;
    iter_preds_test = iter_preds_test+alpha_preds_test;
    
    for (int i = 0; i < y.n_cols; i++) {
      preds_mat.slice(iter).col(i) = iter_preds.col(i) * col_stdev(i);
      preds_mat_test.slice(iter).col(i) = iter_preds_test.col(i) * col_stdev(i);
    }
    for (int i = 0; i < y.n_cols; i++) {
      preds_mat.slice(iter).col(i) += col_means(i);
      preds_mat_test.slice(iter).col(i) += col_means(i);
    }
    for (int i = 0; i < y.n_cols; i++) {
      preds_mat_tau.slice(iter).col(i) = iter_preds_tau.col(i) * col_stdev(i);
      preds_mat_tau_test.slice(iter).col(i) = iter_preds_tau_test.col(i) * col_stdev(i);
    }
    
    sigma=sample_sigma(n, v_0, y_scaled, iter_preds+Z%iter_preds_tau, sigma_0);
    
    arma::mat sigma_scaled = sample_sigma(n, v_0, y, preds_mat.slice(iter)+Z%preds_mat_tau.slice(iter), sigma_0);
    
    sigmas.slice(iter) = sigma_scaled;
    
    ts.slice(iter) = t;
    alpha_means.col(iter) = alpha_mean;
    random_effects.slice(iter) = alpha_is;
    
  }
  
  
  Rcpp::Rcout << "";
  
  
  return List::create(
    Named("predictions") = preds_mat,
    Named("predictions_tau") = preds_mat_tau,
    Named("sigmas") = sigmas,
    Named("ts") = ts,
    Named("alpha_means") = alpha_means,
    Named("random_effects") = random_effects,
    Named("predictions_test") = preds_mat_test,
    Named("predictions_tau_test") = preds_mat_tau_test
  );
}
 

