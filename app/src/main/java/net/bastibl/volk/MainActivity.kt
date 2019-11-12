package net.bastibl.volk

import android.Manifest
import androidx.appcompat.app.AppCompatActivity
import android.os.Bundle
import kotlinx.android.synthetic.main.activity_main.*
import android.content.pm.PackageManager
import android.os.Environment
import androidx.core.app.ActivityCompat
import androidx.core.content.ContextCompat
import java.io.File
import kotlin.concurrent.thread


class MainActivity : AppCompatActivity() {

    private var haveStoragePermission : Boolean = false

    private fun checkStoragePermission() {
        if (ContextCompat.checkSelfPermission(this,
                Manifest.permission.WRITE_EXTERNAL_STORAGE)
            != PackageManager.PERMISSION_GRANTED) {

            ActivityCompat.requestPermissions(this,
                arrayOf(Manifest.permission.WRITE_EXTERNAL_STORAGE),
                123)
        } else {
            haveStoragePermission = true

        }
    }

    override fun onRequestPermissionsResult(requestCode: Int,
                                            permissions: Array<String>, grantResults: IntArray) {
        when (requestCode) {
            123 -> {
                haveStoragePermission = true
            }
        }
    }

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_main)

        checkStoragePermission()

        button.setOnClickListener {
            if (!haveStoragePermission) {
                sample_text.text = "no storage permssion"
                return@setOnClickListener
            }

            button.isEnabled = false;

            thread(start = true) {

                val f = File(Environment.getExternalStorageDirectory(), "/volk/");
                f.mkdirs();

                runOnUiThread {
                    sample_text.text = "benchmarking..."
                }

                val s = runVolk()

                runOnUiThread {
                    sample_text.text = s
                    button.isEnabled = true;
                }
            }
        }
    }

    private external fun runVolk(): String

    companion object {
        init {
            System.loadLibrary("native-lib")
        }
    }
}
