ArrayList<PVector> floodFill(boolean[][] grid, boolean[][] visited, int startX, int startY) {
  ArrayList<PVector> region = new ArrayList<>();
  ArrayList<PVector> stack = new ArrayList<>();
  stack.add(new PVector(startX, startY));

  while (stack.size() > 0) {
    PVector p = stack.remove(stack.size()-1);
    int x = (int)p.x;
    int y = (int)p.y;

    if (x < 0 || y < 0 || x >= grid.length || y >= grid[0].length) continue;
    if (visited[x][y] || grid[x][y]) continue;

    visited[x][y] = true;
    region.add(p);

    stack.add(new PVector(x+1,y));
    stack.add(new PVector(x-1,y));
    stack.add(new PVector(x,y+1));
    stack.add(new PVector(x,y-1));
  }

  return region;
}

boolean[][] keepLargestOpenRegion(boolean[][] grid, int x, int y) {
  boolean[][] visited = new boolean[x][y];
  
  ArrayList<ArrayList<PVector>> regions = new ArrayList<>();

  for(int i=0; i<x ;i++) {
    for(int j=0; j<y ;j++) {
      if (!grid[i][j] && !visited[i][j]) {
        ArrayList<PVector> region = floodFill(grid, visited, i, j);
        regions.add(region);
      }
    }
  }

  ArrayList<PVector> largest = regions.get(0);
  for (ArrayList<PVector> r : regions)
    if (r.size() > largest.size())
      largest = r;


  for (int i=0; i<x ;i++)
    for (int j=0; j<y ;j++)
      grid[i][j] = true;

  for (PVector p : largest)
    grid[(int)p.x][(int)p.y] = false;

  return grid;
}

boolean[][] smoothCave(boolean[][] grid, int x, int y) {
  boolean[][] newGrid = new boolean[x][y];

  for(int i=1; i<x-1 ;i++) {
    for(int j=1; j<y-1 ;j++) {
      int solidCount = 0;
      
      for (int di=-1; di<=1; di++) {
        for (int dj=-1; dj<=1; dj++) {
          if (di == 0 && dj == 0) continue;
          if (grid[i+di][j+dj]) solidCount++;
        }
      }

           if (solidCount > 4) newGrid[i][j] = true;
      else if (solidCount < 4) newGrid[i][j] = false;
      else                     newGrid[i][j] = grid[i][j];
    }
  }

  return newGrid;
}

boolean[][] generateCaves(int x, int y, float fillProbability, int steps) {
  boolean[][] grid = new boolean[x][y];
  
  for(int i=0; i<x ;i++) {
    for(int j=0; j<y ;j++) {
      if (i==0 || j==0 || i==x-1 || j==y-1) grid[i][j] = true;
      else  grid[i][j] = random(1) < fillProbability;
    }
  }

  for(int k=0; k<steps ;k++) grid = smoothCave(grid, x, y);
  
  grid = keepLargestOpenRegion(grid, x, y);
  
  return grid;
}
