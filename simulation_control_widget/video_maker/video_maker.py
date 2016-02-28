"""

created by: lskrinjar
date of creation: 25/02/2016
time of creation: 11:28
"""
import os
import shutil
import threading


try:
    from moviepy.editor import ImageSequenceClip
except:
    ImageSequenceClip = None


from MBD_system.check_filename import check_filename


class VideoMaker(threading.Thread):
    """
    classdocs
    """

    def __init__(self, parent=None):
        """

        :return:
        """
        threading.Thread.__init__(self)

        #   parent
        self._parent = parent
        
        #    name
        self._name = "_animation"


    def run(self):
        """

        :return:
        """
        #   check if folder exists and if not create it
        folder = self._name
        if not os.path.exists(folder):
            os.makedirs(folder)

        #   absolute path to store animation images to use them to create video file
        animation_folder_abs_path = os.path.join(os.getcwd(), folder)

        #   looop through solution data and create figure with opengl function grabFrameBuffer() every i-th step
        for _step in xrange(0, len(self._parent.step), int(self._parent._delta_step)):
            filename = "step_%06d"%_step+".png"

            #   assign step and repaint GL widget
            self._parent._step = int(_step)
            self._parent._update_GL()

            #   get image of current opengl widget scene
            image = self._parent.OpenGLWidget.grabFrameBuffer()

            #   abs path to image object of current simulation time step
            file_abs_path = os.path.join(animation_folder_abs_path, filename)
            #   save image to file
            image.save(file_abs_path)


        #   create video object
        video = ImageSequenceClip(animation_folder_abs_path, fps=24)#24
        #   write data to video object
        __animation_filename = self._name+".avi"

        #   check filename
        __animation_filename = check_filename(__animation_filename)

        video.write_videofile(__animation_filename,codec='mpeg4')

        #   delete animation folder with animation figures
        shutil.rmtree(animation_folder_abs_path)