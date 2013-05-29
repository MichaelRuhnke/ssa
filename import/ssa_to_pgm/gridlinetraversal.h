#ifndef GRIDLINETRAVERSAL_H
#define GRIDLINETRAVERSAL_H
// Sparse Surface Optimization
// Copyright (C) 2011 M. Ruhnke, R. Kuemmerle, G. Grisetti, W. Burgard
//
// SSA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SSA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include <cstdlib>
#include <Eigen/Core>

namespace ssa {

#define GRIDTRAVERSAL_MAXPOINTS 65536

struct GridLineTraversalLine {
  int     num_points;
  Eigen::Vector2i  points[GRIDTRAVERSAL_MAXPOINTS];
};

struct GridLineTraversal {
  inline static void gridLine( Eigen::Vector2i start, Eigen::Vector2i end, GridLineTraversalLine *line ) ;
  inline static void gridLineCore( Eigen::Vector2i start, Eigen::Vector2i end, GridLineTraversalLine *line ) ;

};

void GridLineTraversal::gridLineCore( Eigen::Vector2i start, Eigen::Vector2i end, GridLineTraversalLine *line )
{
  int dx, dy, incr1, incr2, d, x, y, xend, yend, xdirflag, ydirflag;
  int cnt = 0;

  dx = std::abs(end(0)-start(0));
  dy = std::abs(end(1)-start(1));

  if (dy <= dx) {
    d = 2*dy - dx; incr1 = 2 * dy; incr2 = 2 * (dy - dx);
    if (start(0) > end(0)) {
      x = end(0); y = end(1);
      ydirflag = (-1);
      xend = start(0);
    } else {
      x = start(0); y = start(1);
      ydirflag = 1;
      xend = end(0);
    }
    line->points[cnt](0)=x;
    line->points[cnt](1)=y;
    cnt++;
    if (((end(1) - start(1)) * ydirflag) > 0) {
      while (x < xend) {
	x++;
	if (d <0) {
	  d+=incr1;
	} else {
	  y++; d+=incr2;
	}
	line->points[cnt](0)=x;
	line->points[cnt](1)=y;
	cnt++;
      }
    } else {
      while (x < xend) {
	x++;
	if (d <0) {
	  d+=incr1;
	} else {
	  y--; d+=incr2;
	}
	line->points[cnt](0)=x;
	line->points[cnt](1)=y;
	cnt++;
      }
    }
  } else {
    d = 2*dx - dy;
    incr1 = 2*dx; incr2 = 2 * (dx - dy);
    if (start(1) > end(1)) {
      y = end(1); x = end(0);
      yend = start(1);
      xdirflag = (-1);
    } else {
      y = start(1); x = start(0);
      yend = end(1);
      xdirflag = 1;
    }
    line->points[cnt](0)=x;
    line->points[cnt](1)=y;
    cnt++;
    if (((end(0) - start(0)) * xdirflag) > 0) {
      while (y < yend) {
	y++;
	if (d <0) {
	  d+=incr1;
	} else {
	  x++; d+=incr2;
	}
	line->points[cnt](0)=x;
	line->points[cnt](1)=y;
	cnt++;
      }
    } else {
      while (y < yend) {
	y++;
	if (d <0) {
	  d+=incr1;
	} else {
	  x--; d+=incr2;
	}
	line->points[cnt](0)=x;
	line->points[cnt](1)=y;
	cnt++;
      }
    }
  }
  line->num_points = cnt;
}

void GridLineTraversal::gridLine( Eigen::Vector2i start, Eigen::Vector2i end, GridLineTraversalLine *line ) {
  int i,j;
  int half;
  Eigen::Vector2i v;
  gridLineCore( start, end, line );
  if ( start(0)!=line->points[0](0) ||
       start(1)!=line->points[0](1) ) {
    half = line->num_points/2;
    for (i=0,j=line->num_points - 1;i<half; i++,j--) {
      v = line->points[i];
      line->points[i] = line->points[j];
      line->points[j] = v;
    }
  }
}

};//namespace ssa
#endif
