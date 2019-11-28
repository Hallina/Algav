package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.awt.Rectangle;

import supportGUI.Circle;
import supportGUI.Line;

public class DefaultTeam {

	// calculDiametre: ArrayList<Point> --> Line
	// renvoie une pair de points de la liste, de distance maximum.
	public Line calculDiametre(ArrayList<Point> points) {
		if (points.size() < 3) {
			return null;
		}

		Point p = points.get(0);
		Point q = points.get(1);

		double dist = 0;
		double tmp;

		Point t1 = p;
		Point t2 = q;

		for (int i = 0; i < points.size(); i++) {
			t1 = points.get(i);
			for (int j = 1; j < points.size(); j++) {
				t2 = points.get(j);
				tmp = Math.sqrt((t1.y - t2.y) * (t1.y - t2.y) + (t1.x - t2.x) * (t1.x - t2.x));
				if (tmp > dist) {
					dist = tmp;
					p = t1;
					q = t2;
				}
			}
		}

		return new Line(p, q);
	}

	// calculDiametreOptimise: ArrayList<Point> --> Line
	// renvoie une pair de points de la liste, de distance maximum.
	public Line calculDiametreOptimise(ArrayList<Point> points) {
		if (points.size() < 3) {
			return null;
		}

		Point p = points.get(1);
		Point q = points.get(2);

		/*******************
		 * PARTIE A ECRIRE *
		 *******************/
		return new Line(p, q);
	}

	// calculCercleMin: ArrayList<Point> --> Circle
	// renvoie un cercle couvrant tout point de la liste, de rayon minimum.
	public Circle calculCercleMin(ArrayList<Point> points) {
		if (points.isEmpty()) {
			return null;
		}

		Point center = points.get(0);
		int radius = 100;

		/*******************
		 * PARTIE A ECRIRE *
		 *******************/
		return new Circle(center, radius);
	}

	// enveloppeConvexe: ArrayList<Point> --> ArrayList<Point>
	// renvoie l'enveloppe convexe de la liste.
	public ArrayList<Point> enveloppeConvexe(ArrayList<Point> points) {

		if (points.size() < 3) {
			return null;
		}

        return points;
    }
	
	@SuppressWarnings("unchecked")
	public ArrayList<Point> tme6CalculEnveloppe(ArrayList<Point> points) {

		if (points.size() < 3) {
			return null;
		}
		
		if (points.size()<4) return points;

        Point ouest = points.get(0);
        Point sud = points.get(0);
        Point est = points.get(0);
        Point nord = points.get(0);
        for (Point p: points){
            if (p.x<ouest.x) ouest=p;
            if (p.y>sud.y) sud=p;
            if (p.x>est.x) est=p;
            if (p.y<nord.y) nord=p;
        }
        ArrayList<Point> result = new ArrayList<Point>();
        result.add(ouest);
        result.add(sud);
        result.add(est);
        result.add(nord);

        ArrayList<Point> rest = (ArrayList<Point>)points.clone();
        for (int i=0;i<rest.size();i++) {
            if (triangleContientPoint(ouest,sud,est,rest.get(i)) ||
                    triangleContientPoint(ouest,est,nord,rest.get(i))) {
                rest.remove(i);
                i--;
                    }
        }

        for (int i=0;i<result.size();i++) {
            Point a = result.get(i);
            Point b = result.get((i+1)%result.size());
            Point ref = result.get((i+2)%result.size());

            double signeRef = crossProduct(a,b,a,ref);
            double maxValue = 0;
            Point maxPoint = a;

            for (Point p: points) {
                double piki = crossProduct(a,b,a,p);
                if (signeRef*piki<0 && Math.abs(piki)>maxValue) {
                    maxValue = Math.abs(piki);
                    maxPoint = p;
                }java.awt.Rectangle
            }
            if (maxValue!=0){
                for (int j=0;j<rest.size();j++) {
                    if (triangleContientPoint(a,b,maxPoint,rest.get(j))){
                        rest.remove(j);
                        j--;
                    }
                }
                result.add(i+1,maxPoint);
                i--;
            }
        }
        return result;
    }

	private double crossProduct(Point p, Point q, Point s, Point t){
        return ((q.x-p.x)*(t.y-s.y)-(q.y-p.y)*(t.x-s.x));
    }
	
	private boolean triangleContientPoint(Point a, Point b, Point c, Point x) {
        double l1 = ((b.y-c.y)*(x.x-c.x)+(c.x-b.x)*(x.y-c.y))/(double)((b.y-c.y)*(a.x-c.x)+(c.x-b.x)*(a.y-c.y));
        double l2 = ((c.y-a.y)*(x.x-c.x)+(a.x-c.x)*(x.y-c.y))/(double)((b.y-c.y)*(a.x-c.x)+(c.x-b.x)*(a.y-c.y));
        double l3 = 1-l1-l2;
        return (0<l1 && l1<1 && 0<l2 && l2<1 && 0<l3 && l3<1);
    }
	
	public  ArrayList<Point> triPixel(ArrayList<Point> points) {
		int xmax = 0;
		for (Point p : points)
			if (p.x > xmax)
				xmax = p.x;

		Point[] ymin = new Point[xmax + 10];
		Point[] ymax = new Point[xmax + 10];

		for (Point p : points) {
			if (ymin[p.x] == null || ymin[p.x].y > p.y)
				ymin[p.x] = p;
			if (ymax[p.x] == null || ymax[p.x].y < p.y)
				ymax[p.x] = p;
		}
		ArrayList<Point> parcours = new ArrayList<Point>();
		for (int i = 0; i < ymin.length; i++)
			if (ymin[i] != null) {
				parcours.add(ymin[i]);
			}

		for (int i = ymax.length - 1; i >= 0; i--)

			if (ymax[i] != null) {
				parcours.add(ymax[i]);
			}
		
		for (int i = 0; i < parcours.size(); i++)
			if (parcours.get(i).equals(parcours.get((i+1)%parcours.size()))) {
				System.out.println(parcours.get(i));
				parcours.remove(i);
			}
		
		return parcours;
	}

	public Rectangle TME7naif (ArrayList<Point> points){
		ArrayList<Point> enveloppe = tme6CalculEnveloppe(points);
		
		Rectangle rec = new Rectangle();
		
		for (int i = 0; i<enveloppe.size() ; i++) {
			Point p = enveloppe.get(i);
			Point q = enveloppe.get((i + 1)%enveloppe.size());
			
			Rectangle bis = new Rectangle();
			bis = calcul(p, q, enveloppe);
		}
		
		return null;
		
	}

	private double distancePoints(Point a, Point b) {
		
		return Math.sqrt( Math.pow((b.x - a.x), 2) + Math.pow((b.y - a.y), 2)); 
	}
	
	private Rectangle calcul(Point p, Point q, ArrayList<Point> enveloppe) {
		Line l = new Line (p, q);
		
		
//		public double getArea() {
//	        return (side1 + side2 + side3) / 2;
//	    }
		
		return null;
	}

}
