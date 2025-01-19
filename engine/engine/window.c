#include "window.h"

#include "common/status.h"
#include "utils/logger.h"

#include "utils/timer.h"

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch (uMsg)
    {
    case WM_NCCREATE:
    {
        // TODO: Find the best way to do this.

        // Recover the window pointer.
        LPCREATESTRUCT lpcs = (LPCREATESTRUCT)lParam;
        Window* window = (Window*)lpcs->lpCreateParams;

        if (!window)
        {
            log_error("Failed to recover window pointer.");
        }

        // Store the pointer safely for future use.
        SetWindowLongPtrA(hwnd, GWLP_USERDATA, (LONG_PTR)window);
        break;
    }
    case WM_DESTROY:
    {
        PostQuitMessage(0);
        return S_OK;
    }
    case WM_EXITSIZEMOVE:
    {
        Window* window = (Window*)GetWindowLongPtrA(hwnd, GWLP_USERDATA);

        // Calculate the new window dimensions.
        RECT rect;
        GetClientRect(window->hwnd, &rect);

        int width = rect.right - rect.left;
        int height = rect.bottom - rect.top;

        // Check the dimensions are different.
        if (width != window->width || height != window->height)
        {
            // Update the dimensions.
            window->width = width;
            window->height = height;
            
            // TODO: The engine updates the canvas size, which sets the window bitmap
            //       dimensions. Feels off. Not sure.

            // Fire the resize callback.
            window->on_resize(window->ctx);
        }

        break;
    }
    case WM_KEYUP:
    {
        Window* window = (Window*)GetWindowLongPtrA(hwnd, GWLP_USERDATA);
        window->on_keyup(window->ctx, wParam);

        break;
    }
    case WM_INPUT:
    {
        // TODO: I'm not sure if this is causing lag because we're doing it multiple
        //       times a frame. The only reason I'm doing this is to improve performance
        //       but if it's not working, there is no point.
        
        // Use raw input for the mouse, so that we don't have to 
        // reset the mouse position every frame as this consistently
        // took about 1ms, wayyyy too long.
        Window* window = (Window*)GetWindowLongPtrA(hwnd, GWLP_USERDATA);

        UINT dwSize = 0;
        
        // Not sure why this has to be called twice, but doesn't work 
        // without it.
        GetRawInputData((HRAWINPUT)lParam, RID_INPUT, NULL, &dwSize, 
            sizeof(RAWINPUTHEADER));

        char buffer[sizeof(RAWINPUT)] = { 0 };

        if (GetRawInputData((HRAWINPUT)lParam, RID_INPUT, buffer, &dwSize, 
            sizeof(RAWINPUTHEADER)) != dwSize)
        {
            log_error("GetRawInputData did not return the correct size.");
        }
            
        RAWINPUT* raw = (RAWINPUT*)buffer;
        if (raw)
        {
            if (RIM_TYPEMOUSE == raw->header.dwType)
            {
                window->mouse_dx += raw->data.mouse.lLastX;
                window->mouse_dy += raw->data.mouse.lLastY;
            }
        }
        break;
    }
    }

    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

Status window_init(Window* window, Canvas* canvas, void* ctx, int width, int height)
{
	log_info("Initialising the window.");
	memset(window, 0, sizeof(Window));

    window->canvas = canvas;
    window->ctx = ctx;
    window->width = width;
    window->height = height;

    DWORD window_style = WS_OVERLAPPEDWINDOW;

    // Use our desired client area (width, height) to calculate the 
    // full size of the window including the titlebar/borders.
    RECT rect = {
        .right = window->width, 
        .bottom = window->height
    };

    if (!AdjustWindowRect(&rect, window_style, 0))
    {
        log_error("Failed to AdjustWindowRect.");
        return STATUS_WIN32_FAILURE;
    }

    // Use the full size of the window, not the client area.
    const int actual_width = rect.right - rect.left;
    const int actual_height = rect.bottom - rect.top;

    HMODULE hinstance = GetModuleHandleA(NULL);

    // Define the window's class.
    WNDCLASS wc = {
        .lpszClassName = range_WINDOW_CLASS,
        .hInstance = hinstance,
        .lpfnWndProc = WindowProc,
        .hIcon = LoadIcon(NULL, IDI_APPLICATION),
        .hCursor = LoadCursor(NULL, IDC_ARROW),
        .style = CS_OWNDC
    };

    if (!RegisterClassA(&wc))
    {
        log_error("Failed to RegisterClassA.");
        return STATUS_WIN32_FAILURE;
    }

    // Create the window
    window->hwnd = CreateWindowExA(
        0,                          // Window styles, TODO: PASS window_style?
        range_WINDOW_CLASS,         // Window class
        range_WINDOW_TITLE,         // Window caption
        window_style,

        // Size and position
        CW_USEDEFAULT,
        CW_USEDEFAULT,
        actual_width,
        actual_height,

        NULL,                   // Parent window    
        NULL,                   // Menu
        hinstance,              // Instance handle
        window                  // Additional application data
    );

    if (NULL == window->hwnd)
    {
        log_error("Failed to CreateWindowExA.");
        return STATUS_WIN32_FAILURE;
    }

    // Get the Device Context, as we are only drawing to it from this class,
    // I believe it is fine to keep this handle and not release it.
    window->hdc = GetDC(window->hwnd);

    // TODO: Think about this. Can I make it faster.
    // Look at: https://www.youtube.com/watch?v=hNKU8Jiza2g&list=PLnuhp3Xd9PYTt6svyQPyRO_AAuMWGxPzU&index=16
    // Handmade Hero 004.

    // Get the new size and create the frame bitmap info.
    window->bitmap.bmiHeader.biSize = sizeof(window->bitmap.bmiHeader);
    window->bitmap.bmiHeader.biWidth = window->canvas->width;
    window->bitmap.bmiHeader.biHeight = -window->canvas->height;
    window->bitmap.bmiHeader.biPlanes = 1;
    window->bitmap.bmiHeader.biBitCount = 32;
    window->bitmap.bmiHeader.biCompression = BI_RGB; // Uncompressed.

    // Subscribe to mouse raw input.
    RAWINPUTDEVICE raw_input_devices[1] = { 0 };
    raw_input_devices[0].usUsagePage = 0x01; // HID_USAGE_PAGE_GENERIC
    raw_input_devices[0].usUsage = 0x02;     // HID_USAGE_GENERIC_MOUSE
    raw_input_devices[0].dwFlags = 0;
    raw_input_devices[0].hwndTarget = window->hwnd;
    
    if (!RegisterRawInputDevices(raw_input_devices, 1, sizeof(RAWINPUTDEVICE)))
    {
        log_error("Failed to register mouse for raw input.\n");
    }

    ShowWindow(window->hwnd, SW_SHOW);

    return STATUS_OK;
}

int window_process_messages()
{
    // TODO: Some heap is corrupted so i keep getting an error here......
    // TODO: The issue comes when we close the window sometimes and sometimes on startup.
    //       Potentially something to do with storing the window pointer.

    // TODO: The stupid error has to be from some array allocation right? Check the new temp shadow ones. 
    // Error triggered by events.

    // Processes all messages and sends them to WindowProc
    MSG msg;
    while (PeekMessageA(&msg, NULL, 0, 0, PM_REMOVE))
    {
        if (msg.message == WM_QUIT)
        {
            log_info("WM_QUIT message received.");
            return FALSE;
        }

        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return TRUE;
}

void window_display(Window* window)
{
    // TODO: Could look into ways of speeding up blitting, however,
    //       not sure it's really possible, unless there is some 
    //       issue with the width/height needing to be multiples 
    //       of something..

    if (window->canvas->width != window->width || window->canvas->height != window->height)
    {
        StretchDIBits(
            window->hdc,
            0, 0,
            window->width, window->height,
            0, 0,
            window->canvas->width, window->canvas->height,
            window->canvas->pixels,
            &window->bitmap,
            DIB_RGB_COLORS,
            SRCCOPY);
    }
    else
    {
        // This gives slightly better performance than StretchDIBits.
        SetDIBitsToDevice(window->hdc,
            0, 0,
            window->width, window->height,
            0, 0,
            0,
            window->height,
            window->canvas->pixels,
            &window->bitmap,
            DIB_RGB_COLORS);
    }
}

void window_destroy(Window* window)
{
    // TODO: Not sure how necessary this is.
}
