package edu.utah.cs.trajectory

import edu.utah.cs.spatial.{LineSegment, Point}

case class TrajMeta(traj_id: Int, seg_id: Int)

case class Trajectory(id: Int, segments: Array[Point]) {
  def distanceFrom(otherTraj: Trajectory): Double = {
    Math.min(Trajectory.hDistance(this, otherTraj), Trajectory.hDistance(otherTraj, this))
  }
}

object Trajectory {
  var EPSILON = 0.0001
  var DELTA = 3

  def RDPCompress(traj: Array[Point], epsilon: Double): Array[Point] = {
    val baseLineSeg = LineSegment(traj.head, traj.last)
    val dmax = traj.map(x => x.minDist(baseLineSeg)).zipWithIndex.maxBy(_._1)
    if (dmax._1 > epsilon) {
      RDPCompress(traj.slice(0, dmax._2 + 1), epsilon) ++ RDPCompress(traj.slice(dmax._2, traj.length), epsilon)
    } else {
      Array(traj.head, traj.last)
    }
  }

  def parseLine(line: String): Trajectory = {
    val splitted = line.split(" ")
    Trajectory(splitted(0).toInt, splitted.iterator.drop(1).map(_.toDouble).grouped(2).map(seq => Point(Array(seq(0), seq(1)))).toArray)
  }

  def hDistance(traj1: Trajectory, traj2: Trajectory): Double = {
    traj1.segments.iterator.take(traj1.segments.length - 1).zip(traj1.segments.iterator.drop(1)).map {
      case (q0, q1) =>
        val qSegment = LineSegment(q0, q1)
        traj2.segments.iterator.take(traj2.segments.length - 1).zip(traj2.segments.iterator.drop(1)).map {
          case (p0, p1) => qSegment.minDist(LineSegment(p0, p1))
        }.min
    }.max
  }

  def distanceFrom(seg_iter: Iterable[Tuple2[Int, LineSegment]],
      traj2: Array[LineSegment]): Double = {
    seg_iter.map { case (_, seg1) =>
        traj2.iterator.map { seg2 =>
            seg1.minDist(seg2)
        }.min
    }.max
  }

  def hausdorffDistance(x: Array[LineSegment], y: Array[LineSegment]): Double = {
    Math.max(x.map(seg_1 => y.map(seg_2 => seg_1.matchDist(seg_2)).min).max,
             y.map(seg_1 => x.map(seg_2 => seg_1.matchDist(seg_2)).min).max)
  }

  def discreteFrechetDistance(x: Array[LineSegment], y: Array[LineSegment]): Double = {
    val n = x.length
    val m = y.length
    val ca: Array[Array[Double]] = Array.fill[Double](n, m)(-1.0)
    var i = 0
    while (i < n) {
      var j = 0
      while (j < m) {
        if (i == 0 && j == 0) ca(i)(j) = x(i).matchDist(y(j))
        else if (i == 0) ca(i)(j) = Math.max(ca(i)(j - 1), x(i).matchDist(y(j)))
        else if (j == 0) ca(i)(j) = Math.max(ca(i - 1)(j), x(i).matchDist(y(j)))
        else ca(i)(j) = Math.max(Math.min(Math.min(ca(i - 1)(j), ca(i)(j - 1)), ca(i - 1)(j - 1)), x(i).matchDist(y(j)))
        j += 1
      }
      i += 1
    }
    ca.last.last
  }

  def dtwDistance(x: Array[LineSegment], y: Array[LineSegment]): Double = {
    val points1 = x.map(item => item.start) :+ x.last.end
    val points2 = y.map(item => item.start) :+ y.last.end

    val cost = Array.fill[Double](points1.length, points2.length)(Double.MaxValue)

    for (i <- points1.indices) {
      for (j <- points2.indices) {
        cost(i)(j) = points1(i).minDist(points2(j))
        val left = if (i > 0) cost(i - 1)(j) else Double.MaxValue
        val up = if (j > 0) cost(i)(j - 1) else Double.MaxValue
        val diag = if (i > 0 && j > 0) cost(i - 1)(j - 1) else Double.MaxValue
        val last = math.min(math.min(left, up), diag)
        if (i > 0 || j > 0) {
          cost(i)(j) += last
        }
      }
    }

    cost.last.last
  }

  private def subCost(shape1: Point, shape2: Point): Double = {
    if (shape1.minDist(shape2) <= EPSILON) {
      0
    } else {
      1
    }
  }

  def EDRDistance(x: Array[LineSegment], y: Array[LineSegment], threshold: Double = Double.MaxValue): Double = {
    val points1 = x.map(item => item.start) :+ x.last.end
    val points2 = y.map(item => item.start) :+ y.last.end

    val cost = Array.ofDim[Double](points1.length, points2.length)

    for (i <- points1.indices) {
      var minDist = Double.MaxValue
      for (j <- points2.indices) {
        if (i > 0 || j > 0) {
          val left = if (i > 0) cost(i - 1)(j) + 1 else Double.MaxValue
          val up = if (j > 0) cost(i)(j - 1) + 1 else Double.MaxValue
          val diag = if (i > 0 && j > 0) cost(i - 1)(j - 1) + subCost(points1(i), points2(j)) else Double.MaxValue
          val last = math.min(math.min(left, up), diag)
          if (last > threshold) {
            cost(i)(j) = Double.MaxValue
          } else {
            cost(i)(j) = last
          }
        } else {
          cost(i)(j) = subCost(points1(i), points2(j))
        }
        minDist = math.min(minDist, cost(i)(j))
      }

      if (minDist > threshold) {
        return threshold
      }
    }

    cost.last.last
  }

  private def subCost(shape1: Point, index1: Int, shape2: Point, index2: Int): Double = {
    if (math.abs(index1 - index2) <= DELTA && shape1.minDist(shape2) <= EPSILON) {
      0
    } else {
      1
    }
  }

  def LCSSDistance(x: Array[LineSegment], y: Array[LineSegment], threshold: Double = Double.MaxValue): Double = {
    val points1 = x.map(item => item.start) :+ x.last.end
    val points2 = y.map(item => item.start) :+ y.last.end

    val cost = Array.ofDim[Double](points1.length, points2.length)

    for (i <- points1.indices) {
      var minDist = Double.MaxValue
      for (j <- points2.indices) {
        if (i > 0 || j > 0) {
          val left = if (i > 0) cost(i - 1)(j) + 1 else Double.MaxValue
          val up = if (j > 0) cost(i)(j - 1) + 1 else Double.MaxValue
          val diag = if (i > 0 && j > 0) {
            cost(i - 1)(j - 1) + subCost(points1(i), i, points2(j), j)
          } else Double.MaxValue
          val last = math.min(math.min(left, up), diag)
          if (last > threshold) {
            cost(i)(j) = Double.MaxValue
          } else {
            cost(i)(j) = last
          }
        } else {
          cost(i)(j) = subCost(points1(i), i, points2(j), j)
        }
        minDist = math.min(minDist, cost(i)(j))
      }

      if (minDist > threshold) {
        return threshold
      }
    }

    cost.last.last
  }

}
