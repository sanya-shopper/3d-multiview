/* Linux (V4L2) camera streamer for the multiview live protocol.
 * OWNERSHIP (parallel build): this file belongs to the v4l2 work item. */
#ifdef __linux__
#error "work item under construction"
#else
#include <stdio.h>
int main(void)
{
    fprintf(stderr, "stream_cam_v4l2 is Linux-only; use stream_cam "
                    "(Swift) on macOS\n");
    return 1;
}
#endif
