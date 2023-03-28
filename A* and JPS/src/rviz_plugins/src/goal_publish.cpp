#include <tf/transform_listener.h>

#include <geometry_msgs/PoseStamped.h>

#include "rviz/display_context.h"
#include "rviz/properties/string_property.h"

#include "goal_tool.h"
using namespace rviz;

int main(int argc, char *argv[])
{
    setlocale(LC_ALL,"");
    ros::init(argc,argv,"goal_pub");

    return 0;
}
