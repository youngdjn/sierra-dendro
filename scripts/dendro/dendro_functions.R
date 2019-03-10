norm = function(x) {
  
  x_lag = lag(x)
  x_lead = lead(x)
  neighbor_mean = ((x_lag+x_lead)/2)
  
  x_prop = (x-neighbor_mean)/neighbor_mean
  
  return(x_prop)
  
}




get_tree_name = function(core_names) {
  
  tree_names = str_replace(core_names,pattern="[aAbBzZ]",replacement="")
  
  new_names = NULL
  for(i in 1:length(tree_names)) {
    tree_name = tree_names[i]
    new_names[i] = ifelse(tree_name %in% trees$former.id,trees[trees$former.id==tree_name,]$tree.id,tree_name)
  }
  
  return(new_names)
}
