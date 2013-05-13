public class explosion {
    public static void main(String[] args) {
	if (args[0].equals("Boid")) {
	    StdAudio.play("boomBoid.wav");
	    try {
		Thread.sleep(500);
	    } catch (InterruptedException ex) {

	    }
	    System.exit(0);
	} else if (args[0].equals("Player")) {
	    StdAudio.play("boomPlayer.wav");
	    try {
		Thread.sleep(2500);
	    } catch (InterruptedException ex) {

	    }
	    System.exit(0);
	} else if (args[0].equals("Enemy")) {
	    StdAudio.play("boomEnemy.wav");
	    try {
		Thread.sleep(4000);
	    } catch (InterruptedException ex) {

	    }
	    System.exit(0);
	}
    }
}