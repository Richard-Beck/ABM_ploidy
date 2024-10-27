package onLatticeCA_jar;
import HAL.Gui.UIGrid;
import HAL.Gui.UILabel;
import HAL.Gui.UIWindow;

import static HAL.Util.RGB;



public class visualize {
    public UIWindow win;
    public OnLatticeGrid model;

    public visualize(OnLatticeGrid model, String title){

        this.model=model;
        win=new UIWindow(title,true);

        win.AddCol(0,model.vis);
        win.AddCol(1,model.resources.currV);

        win.RunGui();

    }

public void close(){
        win.Close();
}

}
